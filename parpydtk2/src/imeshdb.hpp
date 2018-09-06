// MIT License
//
// Copyright (c) 2018 Qiao Chen
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// \file src/imeshdb.hpp
/// \brief interface mesh database manager

#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <vector>

#include "field.hpp"

// trilinos
#include <DTK_MoabManager.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_MultiVector.hpp>

// moab
#include <MBParallelConventions.h>
#include <moab/ParallelComm.hpp>  // mpi.h in here

namespace parpydtk2 {

/** \addtogroup mesh
 * @{
 */ 

/// \class IMeshDB
/// \brief interface mesh database, build on top of MOAB
class IMeshDB {
  /// \brief helper for clean up mesh
  /// \param[in] del whether or not delete mesh
  inline void init_(bool del = false) {
    if (del) {
      mdb_.delete_mesh();
      locals_.clear();
    }
    ::moab::ErrorCode ret;
    ret = mdb_.tag_get_handle(PARALLEL_GID_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
                              gidtag_);
    handle_moab_error(ret);
    // get partition tag also
    int dv = -1;
    ret = mdb_.tag_get_handle(
        PARALLEL_PARTITION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, parttag_,
        ::moab::MB_TAG_CREAT | ::moab::MB_TAG_SPARSE, &dv);
    handle_moab_error(ret);
  }

  /// \brief handle all vectors
  inline void reset_vecs_() noexcept {
    vsets_.resize(1, root_set);
    locals_.resize(1);
    bboxes_.resize(1);
    gbboxes_.resize(1);
    for (int i = 0; i < 3; ++i)
      bboxes_[0][i] = std::numeric_limits<double>::max();
    for (int i = 3; i < 6; ++i)
      bboxes_[0][i] = std::numeric_limits<double>::min();
  }

  /// \brief initialize empty bounding boxes
  /// \param[in] i index if < 0 then init all
  inline void init_bbox_(int i = -1) noexcept {
    static constexpr double min_ = std::numeric_limits<double>::min();
    static constexpr double max_ = std::numeric_limits<double>::max();
    if (i < 0)
      for (auto &box : bboxes_) {
        int j = 0;
        for (; j < 3; ++j) box[j] = max_;
        for (; j < 6; ++j) box[j] = min_;
      }
    else {
      int j = 0;
      for (; j < 3; ++j) bboxes_[i][j] = max_;
      for (; j < 6; ++j) bboxes_[i][j] = min_;
    }
  }

  /// \brief compute all bounding box
  inline void cmpt_bboxes_() {
    ::moab::ErrorCode ret;
    for (int i = 0; i < (int)vsets_.size(); ++i) {
      const auto &vset = vsets_[i];
      const int sz = size();
      const double *X[3];
      auto &box = bboxes_[i];
      ret = mdb_.get_coords(vset, X[0], X[1], X[2]);
      int j = 0, k = 3, dir = 0;
      if (ret == ::moab::MB_SUCCESS) {
        // opt version
        for (; dir < 3; ++dir) {
          auto minmax = std::minmax_element(X[dir], X[dir] + sz);
          box[j++] = *minmax.first;
          box[k++] = *minmax.second;
        }
      } else {
        // loop 1by1
        double X[3];
        const ::moab::Range &vs = locals_[i];
        for (const auto &v : vs) {
          ret = mdb_.get_coords(&v, 1, X);
          handle_moab_error(ret);
          for (dir = j = 0, k = 3; dir < 3; ++dir, ++j, ++k) {
            box[j] = std::min(box[j], X[dir]);
            box[k] = std::max(box[k], X[dir]);
          }
        }
      }
    }
  }

 public:
  /// \brief constructor with communicator
  /// \param[in] comm communicator
  explicit IMeshDB(MPI_Comm comm = MPI_COMM_WORLD)
      : mdb_(),
        par_(new ::moab::ParallelComm(&mdb_, comm), true),
        vsets_(1, root_set),
        locals_(1),
        bboxes_(1),
        gbboxes_(1),
        created_(false),
        usergid_(false),
        empty_(false),
        has_empty_(false) {
    init_();
    init_bbox_(0);
  }

  virtual ~IMeshDB() = default;

  /// \brief get total ranks
  inline int ranks() const noexcept { return par_->size(); }

  /// \brief get my rank
  inline int rank() const noexcept { return par_->rank(); }

  /// \brief get the communicator
  inline MPI_Comm comm() const noexcept { return par_->comm(); }

  /// \brief get mesh
  ::Teuchos::RCP<moab::ParallelComm> pcomm() const noexcept { return par_; }

  /// \name imesh_py_interface
  ///@{

  /// \brief begin to create mesh
  inline void begin_create() {
    if (created_) {
      show_warning("current mesh will be removed");
      init_(true);
      reset_vecs_();
    }
    created_ = false;
  }

  /// \brief create a new vertex set
  /// \deprecated No longer supports sets
  inline void create_vset() {
    entity_t vset;
    ::moab::ErrorCode ret = mdb_.create_meshset(moab::MESHSET_SET, vset);
    handle_moab_error(ret);
    vsets_.emplace_back(vset);
    locals_.push_back(moab::Range());
    bboxes_.push_back(std::array<double, 6>());
    gbboxes_.push_back(std::array<double, 6>());
    init_bbox_(bboxes_.size() - 1);
  }

  /// \brief create vertices
  /// \param[in] nv number of vertices
  /// \param[in] coords coords values
  inline void create_vertices(int nv, const double *coords) {
    ::moab::Range verts;
    ::moab::ErrorCode ret = mdb_.create_vertices(coords, nv, verts);
    handle_moab_error(ret);
    locals_[0].merge(verts);
  }

  /// \brief extract assigned coordinates
  /// \param[out] coords coordinates
  /// \note coords must be at least n*3 where n is the size of the mesh
  inline void extract_vertices(double *coords) const {
    ::moab::ErrorCode ret = mdb_.get_coords(locals_[0], coords);
    handle_moab_error(ret);
  }

  /// \brief assign global IDs
  /// \param[in] nv number of local vertices
  /// \param[in] gids global IDs
  /// \note \a gids should be one-based indices
  inline void assign_gids(int nv, const int *gids) {
    show_warning_if((entity_t)nv != locals_[0].size(),
                    "global ID size and vertex size do not match");
    moab::ErrorCode ret = mdb_.tag_set_data(gidtag_, locals_[0], gids);
    handle_moab_error(ret);
    usergid_ = true;
  }

  /// \brief extract global IDs
  /// \param[out] gids global IDs
  inline void extract_gids(int *gids) const {
    moab::ErrorCode ret =
        mdb_.tag_get_data(gidtag_, locals_[0], reinterpret_cast<void *>(gids));
    handle_moab_error(ret);
  }

  /// \brief finish manupilating the mesh
  /// \param[in] trivial_gid \a true if we let MOAB to compute the GID
  ///
  /// NOTE that if your input mesh is element-based partition and the vertices
  /// you create are mesh nodes, then you have to specify the correct global
  /// IDs in parallel. However, if the coordinates are face centres, then the
  /// global IDs can be trivially computed by MOAB since there are no shared
  /// entities cross different processes.
  inline void finish_create(bool trivial_gid = true) {
    if (created_) {
      show_warning("a mesh has already been created, do nothing");
      return;
    }
    created_ = true;
    if (trivial_gid) {
      show_warning_if(
          usergid_,
          "the trivial_gid is on in finish_create, your GIDs will be ignored");
      for (const auto &vset : vsets_) {
        moab::ErrorCode ret = par_->assign_global_ids(vset, 0);
        handle_moab_error(ret);
      }
    } else {
      show_warning_if(!usergid_ && (ranks() > 1),
                      "critical!! Global IDs are missing!");
    }
    if (ranks() > 1) {
      // this probably does nothing
      for (const auto &vset : vsets_) {
        ::moab::ErrorCode ret = par_->resolve_shared_ents(vset, 0, 0);
        handle_moab_error(ret);
      }
    }

    // NOTE the following part is to resolve the empty partition issues
    // with DTK by assigning a duplication node from mater process to
    // processes that are empty
    if (rank() == 0)
      for (const auto &vset : vsets_)
        throw_error_if(locals_[vset].size() == 0,
                       "master process cannot have empty partition!");
    for (const auto &vset : vsets_) {
      const entity_t first = locals_[vset][0];
      double dup_coord_and_gid[4];  // encode together
      if (rank() == 0) {
        ::moab::ErrorCode ret = mdb_.get_coords(&first, 1, dup_coord_and_gid);
        handle_moab_error(ret);
        // get global id
        int dup_gid;
        ret = mdb_.tag_get_data(gidtag_, &first, 1, &dup_gid);
        handle_moab_error(ret);
        dup_coord_and_gid[3] = dup_gid;
      }
      if (ranks() > 1) {
        // bcast
        int ret = MPI_Bcast(dup_coord_and_gid, 4, MPI_DOUBLE, 0, par_->comm());
        if (ret != MPI_SUCCESS) {
          // NOTE that the default handler in mpi will directly abort, as a
          // mapper we will not modify the mpi handler (this should be the task
          // of application codes)
          char msg[MPI_MAX_ERROR_STRING];
          int dummy, err_cls;
          MPI_Error_string(ret, msg, &dummy);
          MPI_Error_class(ret, &err_cls);
          std::cerr << "FATAL ERROR! MPI failed with code: " << ret
                    << ", error_class: " << err_cls << ", msg: " << msg
                    << ", rank: " << rank() << ", in " << __FUNCTION__ << " at "
                    << __FILE__ << "():" << __LINE__ << '\n';
          MPI_Abort(MPI_COMM_WORLD, ret);
        }
        // check if "my" mesh is empty
        if (locals_[vset].size() == 0) {
          create_vertices(1, dup_coord_and_gid);
          const int dup_gid = dup_coord_and_gid[3];
          assign_gids(1, &dup_gid);
          empty_ = true;
        } else
          empty_ = false;

        // now notice all procs about empty partition information
        std::vector<int> empty_flags(ranks());
        int my_status = empty_;
        ret = MPI_Allgather(&my_status, 1, MPI_INT, empty_flags.data(), 1,
                            MPI_INT, comm());
        if (ret != MPI_SUCCESS) {
          // NOTE that the default handler in mpi will directly abort, as a
          // mapper we will not modify the mpi handler (this should be the task
          // of application codes)
          char msg[MPI_MAX_ERROR_STRING];
          int dummy, err_cls;
          MPI_Error_string(ret, msg, &dummy);
          MPI_Error_class(ret, &err_cls);
          std::cerr << "FATAL ERROR! MPI failed with code: " << ret
                    << ", error_class: " << err_cls << ", msg: " << msg
                    << ", rank: " << rank() << ", in " << __FUNCTION__ << " at "
                    << __FILE__ << "():" << __LINE__ << '\n';
          MPI_Abort(MPI_COMM_WORLD, ret);
        }
        int counts = std::count(empty_flags.cbegin(), empty_flags.cend(), 1);
        has_empty_ = counts > 0;
        if (has_empty_ && rank() == 0) {
          m2s_.reserve(counts);
          for (int i = 0, N = (int)empty_flags.size(); i < N; ++i)
            if (empty_flags[i]) m2s_.push_back(i);
        }

      } else
        empty_ = false;
    }

    // assign partition based on processes, is this needed in DTK?
    for (const auto &vset : vsets_) {
      int rk = rank();
      ::moab::ErrorCode ret = mdb_.tag_set_data(parttag_, &vset, 1, &rk);
      handle_moab_error(ret);
    }

    // assign some values ot fields
    std::vector<double> values;
    int dim = 1;
    for (const auto &field : fields_) dim = std::max(dim, field.second->dim());
    for (unsigned set_id = 0u; set_id < vsets_.size(); ++set_id) {
      values.resize(size() * dim, 0.0);
      for (auto &field : fields_)
        field.second->assign(locals_[set_id], values.data());
    }

    // set up moab manager
    for (unsigned set_id = 0u; set_id < vsets_.size(); ++set_id) {
      mngrs_.emplace_back(par_, vsets_[set_id], false);
    }
    for (const auto &field : fields_) {
      const std::string &name = *field.second;
      const ::moab::Tag &tag = field.second->tag();
      int set_id = field.second->set();
      dtkfields_.emplace(std::make_pair(
          name, std::make_pair(set_id, mngrs_[set_id].createFieldMultiVector(
                                           vsets_[set_id], tag))));
    }

    // compute bboxes
    cmpt_bboxes_();

    // compute global bounding boxes
    if (ranks() == 1) {
      for (int i = 0; i < (int)bboxes_.size(); ++i)
        for (int j = 0; j < 6; ++j) gbboxes_[i][j] = bboxes_[i][j];
    } else {
      // copy all local first
      for (int i = 0; i < (int)bboxes_.size(); ++i) gbboxes_[i] = bboxes_[i];
      // communication needed
      const int box_len = bboxes_.size() * 6;
      std::vector<double> buffer(box_len * ranks());
      double *sendbuffer;
      bool can_free = false;
      if (bboxes_.size() > 1) {
        sendbuffer = new double[box_len];
        double *pos = sendbuffer;
        for (const auto &box : bboxes_)
          pos = std::copy(box.cbegin(), box.cend(), pos);
        can_free = true;
      } else
        sendbuffer = bboxes_.front().data();
      // allgather
      int ret = MPI_Allgather(sendbuffer, box_len, MPI_DOUBLE, buffer.data(),
                              box_len, MPI_DOUBLE, par_->comm());
      if (can_free) delete[] sendbuffer;
      if (ret != MPI_SUCCESS) {
        // NOTE that the default handler in mpi will directly abort, as a
        // mapper we will not modify the mpi handler (this should be the task
        // of application codes)
        char msg[MPI_MAX_ERROR_STRING];
        int dummy, err_cls;
        MPI_Error_string(ret, msg, &dummy);
        MPI_Error_class(ret, &err_cls);
        std::cerr << "FATAL ERROR! MPI failed with code: " << ret
                  << ", error_class: " << err_cls << ", msg: " << msg
                  << ", rank: " << rank() << ", in " << __FUNCTION__ << " at "
                  << __FILE__ << "():" << __LINE__ << '\n';
        MPI_Abort(MPI_COMM_WORLD, ret);
      }

      // compute the global bounding boxes
      for (int i = 0; i < (int)gbboxes_.size(); ++i) {
        auto &gbox = gbboxes_[i];
        const int ld = i * 6;
        for (int j = 0; j < ranks(); ++j)
          for (int k = 0; k < 3; ++k) {
            gbox[k] = std::min(gbox[k], buffer[box_len * j + ld + k]);
            gbox[k + 3] =
                std::max(gbox[k + 3], buffer[box_len * j + ld + k + 3]);
          }
      }
    }
  }

  /// \brief check if empty partition
  inline bool empty() const noexcept { return empty_; }

  /// \brief check if any of the process has an empty partition
  inline bool has_empty() const noexcept { return has_empty_; }

  /// \brief get a reference to the m2s pattern
  /// \note This is used in Python level as "private" thus having "_"
  inline const std::vector<int> &_m2s() const noexcept { return m2s_; }

  /// \brief check mesh size
  inline int size() const noexcept { return locals_[0].size(); }

  /// \brief get bounding box
  /// \param[out] v values
  inline void get_bbox(double *v) const noexcept {
    const auto &box = bboxes_[0];
    for (int i = 0; i < 6; ++i) v[i] = box[i];
  }

  /// \brief get global bounding box
  /// \param[out] v values
  inline void get_gbbox(double *v) const noexcept {
    const auto &box = gbboxes_[0];
    for (int i = 0; i < 6; ++i) v[i] = box[i];
  }

  /// \brief create a field
  /// \param[in] field_name field name
  /// \param[in] dim field dimension
  inline void create_field(const std::string &field_name, int dim = 1) {
    throw_error_if(dim < 1, "invalid field dimension");
    throw_error_if(
        created_,
        "you cannot create fields once an IMeshDB is marked as created!");
    fields_.create(mdb_, field_name, dim);
  }

  /// \brief check if we have a field
  /// \param[in] field_name field name
  inline bool has_field(const std::string &field_name) const noexcept {
    return fields_.has_field(field_name);
  }

  /// \brief check field dimension
  /// \param[in] field_name field name
  inline int field_dim(const std::string &field_name) const {
    return fields_[field_name].dim();  // operator[] throws
  }

  /// \brief assign a value to a field
  /// \param[in] field_name field name
  /// \param[in] values field data values
  inline void assign_field(const std::string &field_name,
                           const double *values) {
    throw_error_if(!created_,
                   "you cannot assign/extract values on an incomplete mesh.. "
                   "did you forget to call finish_create?");
    fields_[field_name].assign(locals_[0], values);  // operator[] throws
  }

  /// \brief assign to the first node
  /// \param[in] field_name field name
  /// \param[in] values field data values, at least size of field dimension
  ///
  /// This function is used by Python to resolve the issues when empty
  /// empty partitions happen. Therefore, this function has an "_" prefix to
  /// indicate "private" usage!
  inline void _assign_1st(const std::string &field_name, const double *values) {
    fields_[field_name].assign_1st(locals_[0], values);  // operator[] throws
  }

  /// \brief extract value
  /// \param[in] field_name field name
  /// \param[out] values field data values
  inline void extract_field(const std::string &field_name,
                            double *values) const {
    throw_error_if(!created_,
                   "you cannot assign/extract values on an incomplete mesh.. "
                   "did you forget to call finish_create?");
    fields_[field_name].extract(locals_[0], values);  // operator[] throws
  }

  /// \brief extract first value
  /// \param[in] field_name field name
  /// \param[out] values field data values
  ///
  /// This function is used by Python to resolve the issues when empty
  /// empty partitions happen. Therefore, this function has an "_" prefix to
  /// indicate "private" usage!
  inline void _extract_1st(const std::string &field_name,
                           double *values) const {
    fields_[field_name].extract_1st(locals_[0], values);  // operator[] throws
  }

  ///@}

  /// \brief get the manger
  inline std::vector<DataTransferKit::MoabManager> &mangers() noexcept {
    return mngrs_;
  }

  /// \brief set geometry dimension
  /// \param[in] dim dimension
  inline void set_dimension(int dim) {
    throw_error_if(dim < 1 || dim > 3, "invalid dimension");
    ::moab::ErrorCode ret = mdb_.set_dimension(dim);
    handle_moab_error(ret);
  }

  /// \brief check if ready
  inline bool ready() const noexcept { return created_; }

 protected:
  /// \brief moab instance
  ::moab::Core mdb_;

  /// \brief moab parallel interface
  ::Teuchos::RCP<moab::ParallelComm> par_;

  /// \brief vertex sets
  std::vector<entity_t> vsets_;

  /// \brief local vertex range
  std::vector<moab::Range> locals_;

  /// \brief bounding boxes
  std::vector<std::array<double, 6>> bboxes_;

  /// \brief global bounding boxes
  std::vector<std::array<double, 6>> gbboxes_;

  /// \brief field data set
  FieldDataSet fields_;

  /// \brief flag to indicate whether users are done with creating mesh
  bool created_;

  /// \brief global ID tag
  ::moab::Tag gidtag_;

  /// \brief partition tag
  ::moab::Tag parttag_;

  /// \brief flag to indicate if we have user computed global ID
  bool usergid_;

  /// \brief DTK MOAB manager
  std::vector<DataTransferKit::MoabManager> mngrs_;

  /// \brief handy typedef
  typedef std::unordered_map<
      std::string,
      std::pair<int, Teuchos::RCP<Tpetra::MultiVector<
                         double, int, DataTransferKit::SupportId>>>>
      dtk_field_t;

  /// \brief DTK multi vector for MOAB tags
  dtk_field_t dtkfields_;

  /// \brief check empty partition
  bool empty_;

  /// \brief check if any process is empty
  bool has_empty_;

  /// \brief comm pattern for master2slaves for handling empty partitions
  std::vector<int> m2s_;

 public:
  /// \brief get the dtk fields
  inline dtk_field_t &dtk_fields() noexcept { return dtkfields_; }
};

/** @}*/

}  // namespace parpydtk2
