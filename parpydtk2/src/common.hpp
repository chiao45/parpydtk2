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

/// \file src/common.hpp
/// \brief common interface
/// \note exceptions need to be handled properly in parallel

#pragma once

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <moab/Core.hpp>

namespace parpydtk2 {

/** \addtogroup common
 * @{
 */

/// \typedef entity_t
/// \brief MOAB entity handle
typedef ::moab::EntityHandle entity_t;

/// \brief root set
enum { root_set = 0 };

/// \def handle_moab_error
/// \brief macro to handle moab error
#define handle_moab_error(__ret)                                       \
  do {                                                                 \
    if (__ret != ::moab::MB_SUCCESS) {                                 \
      std::ostringstream msg;                                          \
      msg << "ERROR! MOAB failed with error code: [" << __ret << "]|[" \
          << ::moab::ErrorCodeStr[__ret] << "], in " << __FUNCTION__   \
          << "() at " << __FILE__ << ':' << __LINE__;                  \
      throw std::runtime_error(msg.str());                             \
    }                                                                  \
  } while (false)

/// \def throw_error
/// \brief throw runtime_error exception
#define throw_error(__msg)                                           \
  do {                                                               \
    std::ostringstream msg;                                          \
    msg << "ERROR! " << __msg << ", in " << __FUNCTION__ << "() at " \
        << __FILE__ << ':' << __LINE__;                              \
    throw std::runtime_error(msg.str());                             \
  } while (false)

/// \def throw_error_if
/// \brief conditionally error throw
#define throw_error_if(__cond, __msg) \
  if ((__cond)) throw_error(__msg)

/// \def throw_noimpl
/// \brief throw not implemented feature error
#define throw_noimpl(__what)                                           \
  do {                                                                 \
    std::ostringstream msg;                                            \
    msg << "Sorry! " << __what << " is not supported yet!"             \
        << " This message was throwed in " << __FUNCTION__ << "() at " \
        << __FILE__ << ':' << __LINE__;                                \
    throw std::runtime_error(msg.str());                               \
  } while (false)

/// \def throw_noimpl_if
/// \brief throw not implemented feature error with condition
#define throw_noimpl_if(__cond, __what) \
  if ((__cond)) throw_noimpl(__what)

/// \def show_warning(__msg)
/// \brief log warning message in stderr
#define show_warning(__msg)                                                  \
  do {                                                                       \
    std::cerr << "WARNING! " << __msg << ", in " << __FUNCTION__ << "() at " \
              << __FILE__ << ':' << __LINE__ << '\n';                        \
                                                                             \
  } while (false)

/// \def show_warning_if
/// \brief log warning with condition
#define show_warning_if(__cond, __msg) \
  if ((__cond)) show_warning(__msg)

/// \def show_experimental
/// \brief log experimental warning in stderr
#define show_experimental(__msg)                                              \
  do {                                                                        \
    std::cerr << "EXPERIMENTAL WARNING! " << __msg << ", in " << __FUNCTION__ \
              << "() at " << __FILE__ << ':' << __LINE__ << '\n';             \
                                                                              \
  } while (false)

/// \def show_experimental_if
/// \brief log experimental warning in stderr with condition
#define show_experimental_if(__cond, __msg) \
  if ((__cond)) show_experimental(__msg)

/// \def show_info
/// \brief show information in parallel
#define show_info(__msg, __rank)                         \
  do {                                                   \
    std::cout << '[' << __rank << "] " << __msg << '\n'; \
  } while (false)

/// \def show_info_master
/// \brief show information only on master rank
#define show_info_master(__msg, __rank) \
  if (!__rank) std::cout << __msg << '\n'

/// \def streamer
/// \brief streaming message with specific rank
#define streamer(__rank) std::cout << '[' << __rank << "] "

/// \def streamer_master
/// \brief streaming messages only on the master process
#define streamer_master(__rank) \
  if (!__rank) std::cout

/** @}*/

}  // namespace parpydtk2
