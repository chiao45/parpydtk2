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

/// \file src/field.hpp
/// \brief wrapper of MOAB tag for field data

#pragma once

#include <memory>
#include <unordered_map>

#include "common.hpp"

namespace parpydtk2 {

/** \addtogroup field
 * @{
 */

/// \class FieldData
/// \brief a representation of MOAB tag for field data
class FieldData {
 public:
  FieldData() = delete;

  /// \brief constructor with moab instance
  /// \param[in] mdb moab instance
  /// \param[in] field_name field name
  /// \param[in] dim field dimension
  /// \note \a set is not the same as mesh set in MOAB
  FieldData(::moab::Core &mdb, const std::string &field_name, int dim = 1)
      : mdb_(mdb), fn_(field_name), set_(0), dim_(dim) {
    ::moab::ErrorCode ret;
    // create tag
    double dv[dim];
    for (int i = 0; i < dim; ++i) dv[i] = 0.0;
    ret = mdb_.tag_get_handle(fn_.c_str(), dim_, moab::MB_TYPE_DOUBLE, tag_,
                              ::moab::MB_TAG_DENSE | ::moab::MB_TAG_CREAT, dv);
    if (ret == ::moab::MB_ALREADY_ALLOCATED) {
      show_warning(fn_ + " has already be created");
    } else {
      handle_moab_error(ret);
    }
  }

  virtual ~FieldData() = default;

  /// \brief assign values
  /// \param[in] range entity ranges, for this work, it should be vertices
  /// \param[in] values data values, for vector/tensor, C order is expected
  inline void assign(const ::moab::Range &range, const double *values) {
    ::moab::ErrorCode ret =
        mdb_.tag_set_data(tag_, range, reinterpret_cast<const void *>(values));
    handle_moab_error(ret);
  }

  /// \brief assign to first node
  ///
  /// This function is used by Python for handling empty partitions
  ///
  /// \param[in] range entity ranges
  /// \param[in] values values for the first node, at least size of dim
  inline void assign_1st(const ::moab::Range &range, const double *values) {
    const entity_t first = range[0];
    ::moab::ErrorCode ret = mdb_.tag_set_data(tag_, &first, 1, values);
    handle_moab_error(ret);
  }

  /// \brief extract values
  /// \param[in] range entity ranges, for this work, it should be vertices
  /// \param[out] values data values, for vector/tensor, C order is expected
  inline void extract(const ::moab::Range &range, double *values) const {
    moab::ErrorCode ret =
        mdb_.tag_get_data(tag_, range, reinterpret_cast<void *>(values));
    handle_moab_error(ret);
  }

  /// \brief extract the value from the first node
  /// \param[in] range entity ranges
  /// \param[out] values values for the first node, at least size of dim
  inline void extract_1st(const ::moab::Range &range, double *values) const {
    const entity_t first = range[0];
    moab::ErrorCode ret = mdb_.tag_get_data(tag_, &first, 1, values);
    handle_moab_error(ret);
  }

  /// brief implicitly cast to string
  inline operator const std::string &() const noexcept { return fn_; }

  /// \brief check the dimension
  inline int dim() const noexcept { return dim_; }

  /// \brief check the set ID
  inline int set() const noexcept { return set_; }

  /// \brief get MOAB tag
  inline const ::moab::Tag &tag() const noexcept { return tag_; }

 protected:
  /// \brief reference to moab instance
  ::moab::Core &mdb_;

  /// \brief field name
  std::string fn_;

  /// \brief set count
  /// \deprecated not used, always 0
  int set_;

  /// \brief field dimension
  int dim_;

  /// \brief moab tag
  ::moab::Tag tag_;
};

/// \class FieldDataSet
/// \brief a set of field data
class FieldDataSet {
  typedef std::unordered_map<std::string, FieldData *> base_t;
  ///< data structure

 public:
  typedef base_t::iterator iterator;              ///< interator type
  typedef base_t::const_iterator const_iterator;  ///< constant iterator

  /// \brief destructor
  virtual ~FieldDataSet() {
    for (auto iter = fs_.begin(); iter != fs_.end(); ++iter)
      delete iter->second;
  }
  /// \brief check if a field exist
  /// \param[in] fn field name
  inline bool has_field(const std::string &fn) const noexcept {
    return fs_.find(fn) != fs_.cend();
  }

  /// \brief create an data field
  /// \param[in] mdb moab data base
  /// \param[in] field_name field name
  /// \param[in] dim field dimension
  /// \note \a set is not the same as mesh set in MOAB
  inline void create(::moab::Core &mdb, const std::string &field_name,
                     int dim = 1) {
    // use move
    fs_.emplace(
        std::make_pair(field_name, new FieldData(mdb, field_name, dim)));
  }

  /// \brief get a reference to a field data
  /// \param[in] fn field name
  /// \note this overloads the base operator[]
  FieldData &operator[](const std::string &fn) {
    try {
      return *fs_.at(fn);
    } catch (...) {
      throw_error(fn + " did not exist");
    }
  }

  /// \brief get a const reference to a field data
  /// \param[in] fn field name
  const FieldData &operator[](const std::string &fn) const {
    try {
      return *fs_.at(fn);
    } catch (...) {
      throw_error(fn + " did not exist");
    }
  }

  /// \brief get the first iterator
  inline iterator begin() noexcept { return fs_.begin(); }

  /// \brief get the end iterator
  inline iterator end() noexcept { return fs_.end(); }

  /// \brief get the constant iterator
  inline const_iterator begin() const noexcept { return fs_.begin(); }

  /// \brief get the constant end iterator
  inline const_iterator end() const noexcept { return fs_.end(); }

  /// \brief get the constant iterator
  inline const_iterator cbegin() const noexcept { return fs_.cbegin(); }

  /// \brief get the constant end iterator
  inline const_iterator cend() const noexcept { return fs_.cend(); }

 protected:
  base_t fs_;  ///< fields
};

/** @}*/

}  // namespace parpydtk2
