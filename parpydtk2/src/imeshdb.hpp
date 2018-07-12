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

#include <string>

#include "field.hpp"

#include <Teuchos_RCP.hpp>
#include <moab/ParallelComm.hpp>
#include <MBParallelConventions.h>

namespace parpydtk2 {

/// \class IMeshDB
/// \brief interface mesh database, build on top of MOAB
class IMeshDB {
  /// \brief helper for clean up mesh
  /// \param[in] del whether or not delete mesh
  void init_(bool del = false) {
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

public:
  /// \brief constructor with communicator
  /// \param[in] name user defined mesh name
  /// \param[in] comm communicator
  explicit IMeshDB(const std::string &name, MPI_Comm comm = MPI_COMM_WORLD)
      : name_(name), mdb_(), par_(new ::moab::ParallelComm(&mdb_, comm), true),
        vset_(root_set), created_(false), usergid_(false) {
    init_();
  }

  virtual ~IMeshDB() = default;

  /// \brief get total ranks
  int ranks() const noexcept { return par_->size(); }

  /// \brief get my rank
  int rank() const noexcept { return par_->rank(); }

  /// \brief get mesh
  ::Teuchos::RCP<moab::ParallelComm> pcomm() const noexcept { return par_; }

  /// \name imesh_py_interface
  ///@{

  /// \brief begin to create mesh
  inline void begin_create() {
    if (created_) {
      show_warning("current mesh will be removed");
      init_(true);
    }
    created_ = false;
  }

  /// \brief create vertices
  /// \param[in] nv number of vertices
  /// \param[in] coords coords values
  inline void create_vertices(int nv, const double *coords) {
    ::moab::Range     verts;
    ::moab::ErrorCode ret = mdb_.create_vertices(coords, nv, verts);
    handle_moab_error(ret);
    locals_.merge(verts);
  }

  /// \brief assign global IDs
  /// \param[in] nv number of local vertices
  /// \param[in] gids global IDs
  /// \note \a gids should be one-based indices
  inline void assign_gids(int nv, const int *gids) {
    show_warning_if(nv != locals_.size(),
                    "global ID size and vertex size do not match");
    moab::ErrorCode ret = mdb_.tag_set_data(gidtag_, locals_, gids);
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
      moab::ErrorCode ret = par_->assign_global_ids(vset_, 0);
      handle_moab_error(ret);
    } else {
      show_warning_if(!usergid_ && (ranks() > 1),
                      "critical!! Global IDs are missing!");
    }
    if (rank() > 1) {
      // this probably does nothing
      ::moab::ErrorCode ret = par_->resolve_shared_ents(vset_, 0, 0);
      handle_moab_error(ret);
    }
  }

  /// \brief check mesh size
  inline int size() const noexcept { return locals_.size(); }

  /// \brief create a field
  /// \param[in] field_name field name
  /// \param[in] dim field dimension
  inline void create_field(const std::string &field_name, int dim = 1) {
    throw_error_if(dim < 1, "invalid field dimension");
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
    return fields_[field_name].dim(); // operator[] throws
  }

  /// \brief assign a value to a field
  /// \param[in] field_name field name
  /// \param[in] values field data values
  inline void assign_field(const std::string &field_name,
                           const double *     values) {
    throw_error_if(!created_,
                   "you cannot assign/extract values on a incomplete mesh.. "
                   "did you forget to call finish_create?");
    fields_[field_name].assign(locals_, values); // operator[] throws
  }

  /// \brief extract value
  /// \param[in] field_name field name
  /// \param[out] values field data values
  inline void extract_field(const std::string &field_name,
                            double *           values) const {
    throw_error_if(!created_,
                   "you cannot assign/extract values on a incomplete mesh.. "
                   "did you forget to call finish_create?");
    fields_[field_name].extract(locals_, values); // operator[] throws
  }

  ///@}

protected:
  /// \brief mesh name
  std::string name_;

  /// \brief moab instance
  ::moab::Core mdb_;

  /// \brief moab parallel interface
  ::Teuchos::RCP<moab::ParallelComm> par_;

  /// \brief vertice set
  entity_t vset_;

  /// \brief local vertex range
  ::moab::Range locals_;

  /// \brief field data set
  FieldDataSet fields_;

  /// \brief flag to indicate whether users are done with creating mesh
  bool created_;

  /// \brief global ID and partition tags
  ::moab::Tag gidtag_, parttag_;

  /// \brief flag to indicate if we have user computed global ID
  bool usergid_;
};

} // namespace parpydtk2
