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
public:
  /// \brief constructor with communicator
  /// \param[in] name user defined mesh name
  /// \param[in] comm communicator
  explicit IMeshDB(const std::string &name, MPI_Comm comm = MPI_COMM_WORLD)
      : name_(name), mdb_(), par_(new ::moab::ParallelComm(&mdb_, comm), true),
        vset_(root_set), created_(false) {
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
  void begin_create() { throw_noimpl_if(created_, "recreating a mesh"); }

  /// \brief create vertices
  /// \param[in] nv number of vertices
  /// \param[in] coords coords values
  void create_vertices(int nv, const double *coords) {
    ::moab::Range     verts;
    ::moab::ErrorCode ret = mdb_.create_vertices(coords, nv, verts);
    handle_moab_error(ret);
    locals_.merge(verts);
  }

  /// \brief assign global IDs
  /// \param[in] nv number of local vertices
  /// \param[in] gids global IDs
  /// \note \a gids should be one-based indices
  void assign_gids(int nv, const int *gids) {
    show_warning_if(nv != locals_.size(),
                    "global ID size and vertex size do not match");
    moab::ErrorCode ret = mdb_.tag_set_data(gidtag_, locals_, gids);
    handle_moab_error(ret);
  }

  /// \brief finish manupilating the mesh
  void finish_create() {
    if (created_) {
      show_warning("a mesh has already been created, do nothing");
      return;
    }
    created_ = true;
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
};

} // namespace parpydtk2
