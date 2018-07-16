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

/// \file parpydtk2/src/imeshdb.hpp
/// \brief interface mesh database manager

#pragma once

#include <array>
#include <vector>
#include <limits>
#include <algorithm>

#include "field.hpp"

// trilinos
#include <Teuchos_RCP.hpp>
#include <DTK_MoabManager.hpp>
#include <Tpetra_MultiVector.hpp>

// moab
#include <moab/ParallelComm.hpp> // mpi.h in here
#include <MBParallelConventions.h>

namespace parpydtk2 {

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
    // get partiiton tag also
    int dv = -1;
    ret    = mdb_.tag_get_handle(
        PARALLEL_PARTITION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, parttag_,
        ::moab::MB_TAG_CREAT | ::moab::MB_TAG_SPARSE, &dv);
    handle_moab_error(ret);
  }

  /// \brief handle all vectors
  inline void reset_vecs_() noexcept {
    vsets_.resize(1, root_set);
    locals_.resize(1);
    bboxes_.resize(1);
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
        for (; j < 3; ++j)
          box[j] = max_;
        for (; j < 6; ++j)
          box[j] = min_;
      }
    else {
      int j = 0;
      for (; j < 3; ++j)
        bboxes_[i][j] = max_;
      for (; j < 6; ++j)
        bboxes_[i][j] = min_;
    }
  }

  /// \brief compute all bounding box
  inline void cmpt_bboxes_() {
    ::moab::ErrorCode ret;
    for (int i = 0; i < vsets_.size(); ++i) {
      const auto &  vset = vsets_[i];
      const int     sz   = size(i);
      const double *X[3];
      auto &        box = bboxes_[i];
      ret               = mdb_.get_coords(vset, X[0], X[1], X[2]);
      int j = 0, k = 3, dir = 0;
      if (ret == ::moab::MB_SUCCESS) {
        // opt version
        for (; dir < 3; ++dir) {
          auto minmax = std::minmax_element(X[dir], X[dir] + sz);
          box[j++]    = *minmax.first;
          box[k++]    = *minmax.second;
        }
      } else {
        // loop 1by1
        double               X[3];
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
      : mdb_(), par_(new ::moab::ParallelComm(&mdb_, comm), true),
        vsets_(1, root_set), locals_(1), bboxes_(1), created_(false),
        usergid_(false) {
    init_();
    init_bbox_(0);
  }

  virtual ~IMeshDB() = default;

  /// \brief get total ranks
  inline int ranks() const noexcept { return par_->size(); }

  /// \brief get my rank
  inline int rank() const noexcept { return par_->rank(); }

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
  inline void create_vset() {
    entity_t          vset;
    ::moab::ErrorCode ret = mdb_.create_meshset(moab::MESHSET_SET, vset);
    handle_moab_error(ret);
    vsets_.emplace_back(vset);
    locals_.push_back(moab::Range());
    bboxes_.push_back(std::array<double, 6>());
    init_bbox_(bboxes_.size() - 1);
  }

  /// \brief create vertices
  /// \param[in] nv number of vertices
  /// \param[in] coords coords values
  /// \param[in] set_id set ID
  inline void create_vertices(int nv, const double *coords,
                              unsigned set_id = 0u) {
    throw_error_if(set_id >= vsets_.size(), "exceed set count");
    ::moab::Range     verts;
    ::moab::ErrorCode ret = mdb_.create_vertices(coords, nv, verts);
    handle_moab_error(ret);
    locals_[set_id].merge(verts);
  }

  /// \brief assign global IDs
  /// \param[in] nv number of local vertices
  /// \param[in] gids global IDs
  /// \param[in] set_id set ID
  /// \note \a gids should be one-based indices
  inline void assign_gids(int nv, const int *gids, unsigned set_id = 0u) {
    throw_error_if(set_id >= vsets_.size(), "exceed set count");
    show_warning_if(nv != locals_[set_id].size(),
                    "global ID size and vertex size do not match");
    moab::ErrorCode ret = mdb_.tag_set_data(gidtag_, locals_[set_id], gids);
    handle_moab_error(ret);
    usergid_ = true;
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
    if (rank() > 1) {
      // this probably does nothing
      for (const auto &vset : vsets_) {
        ::moab::ErrorCode ret = par_->resolve_shared_ents(vset, 0, 0);
        handle_moab_error(ret);
      }
    }

    // assign partition based on processes, is this needed in DTK?
    for (const auto &vset : vsets_) {
      int               rk  = rank();
      ::moab::ErrorCode ret = mdb_.tag_set_data(parttag_, &vset, 1, &rk);
      handle_moab_error(ret);
    }

    // assign some values ot fields
    std::vector<double> values;
    int                 dim = 1;
    for (const auto &field : fields_)
      dim = std::max(dim, field.second->dim());
    for (unsigned set_id = 0u; set_id < vsets_.size(); ++set_id) {
      values.resize(size(set_id) * dim, 0.0);
      for (auto &field : fields_)
        field.second->assign(locals_[set_id], values.data());
    }

    // set up moab manager
    for (unsigned set_id = 0u; set_id < vsets_.size(); ++set_id) {
      mngrs_.emplace_back(par_, vsets_[set_id], false);
    }
    for (const auto &field : fields_) {
      const std::string &name   = *field.second;
      const ::moab::Tag &tag    = field.second->tag();
      int                set_id = field.second->set();
      dtkfields_.emplace(std::make_pair(
          name, std::make_pair(set_id, mngrs_[set_id].createFieldMultiVector(
                                           vsets_[set_id], tag))));
    }

    // compute bboxes
    cmpt_bboxes_();
  }

  /// \brief check mesh size
  /// \param[in] set_id set ID
  inline int size(unsigned set_id = 0u) const {
    throw_error_if(set_id >= vsets_.size(), "exceed set count");
    return locals_[set_id].size();
  }

  /// \brief check the set number
  inline int sets() const noexcept { return vsets_.size(); }

  /// \brief get bounding box
  /// \param[out] v values
  /// \param[in] set_id set ID
  inline void get_bbox(double *v, unsigned set_id = 0u) const {
    throw_error_if(set_id >= vsets_.size(), "exceed set count");
    const auto &box = bboxes_[set_id];
    for (int i = 0; i < 6; ++i)
      v[i] = box[i];
  }

  /// \brief create a field
  /// \param[in] field_name field name
  /// \param[in] set_id set ID
  /// \param[in] dim field dimension
  inline void create_field(const std::string &field_name, unsigned set_id = 0u,
                           int dim = 1) {
    throw_error_if(dim < 1, "invalid field dimension");
    throw_error_if(set_id >= vsets_.size(), "exceed set count");
    fields_.create(mdb_, field_name, set_id, dim);
  }

  /// \brief check if we have a field
  /// \param[in] field_name field name
  inline bool has_field(const std::string &field_name) const noexcept {
    return fields_.has_field(field_name);
  }

  /// \brief check field dimension
  /// \param[in] field_name field name
  inline int field_dim(const std::string &field_name) const {
    return fields_[field_name].dim(); // operator[] throws
  }

  /// \brief check field set id
  /// \param[in] field_name field name
  inline int field_set_id(const std::string &field_name) const {
    return fields_[field_name].set(); // operator[] throws
  }

  /// \brief assign a value to a field
  /// \param[in] field_name field name
  /// \param[in] values field data values
  /// \param[in] set_id set ID
  inline void assign_field(const std::string &field_name, const double *values,
                           unsigned set_id = 0u) {
    throw_error_if(!created_,
                   "you cannot assign/extract values on a incomplete mesh.. "
                   "did you forget to call finish_create?");
    fields_[field_name].assign(locals_[set_id], values); // operator[] throws
  }

  /// \brief extract value
  /// \param[in] field_name field name
  /// \param[out] values field data values
  /// \param[in] set_id set ID
  inline void extract_field(const std::string &field_name, double *values,
                            unsigned set_id = 0u) const {
    throw_error_if(!created_,
                   "you cannot assign/extract values on a incomplete mesh.. "
                   "did you forget to call finish_create?");
    fields_[field_name].extract(locals_[set_id], values); // operator[] throws
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
  std::vector<std::array<double, 6> > bboxes_;

  /// \brief field data set
  FieldDataSet fields_;

  /// \brief flag to indicate whether users are done with creating mesh
  bool created_;

  /// \brief global ID and partition tags
  ::moab::Tag gidtag_, parttag_;

  /// \brief flag to indicate if we have user computed global ID
  bool usergid_;

  /// \brief DTK MOAB manager
  std::vector<DataTransferKit::MoabManager> mngrs_;

  /// \brief handy typedef
  typedef std::unordered_map<
      std::string,
      std::pair<int, Teuchos::RCP<Tpetra::MultiVector<
                         double, int, DataTransferKit::SupportId> > > >
      dtk_field_t;

  /// \brief DTK multi vector for MOAB tags
  dtk_field_t dtkfields_;

public:
  /// \brief get the dtk fields
  inline dtk_field_t &dtk_fields() noexcept { return dtkfields_; }
};

} // namespace parpydtk2
