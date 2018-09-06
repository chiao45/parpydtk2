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

/// \file src/dtk2.hpp
/// \brief mapper based on dtk2

#pragma once

#include <algorithm>
#include <map>
#include <memory>
#include <set>

#include "imeshdb.hpp"

// trilinos
#include <DTK_FunctionSpace.hpp>
#include <DTK_MapOperatorFactory.hpp>
#include <Teuchos_ParameterList.hpp>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// helper macro
#define FOR_DIR(__i) for (int __i = 0; i < 2; ++i)

#if !defined(LEN1) || (LEN1 <= 0)
// length of '-----' string
#define LEN1 100
#endif
#if !defined(LEN2) || (LEN2 <= 0)
//  white spaces before option names
#define LEN2 15
#endif

#endif

namespace parpydtk2 {

/** \addtogroup mapper
 * @{
 */ 

/// \enum Methods
/// \brief available methods
enum Methods {
  MMLS = 0,  ///< modified moving least square, default
  SPLINE,    ///< spline interpolation
  N2N        ///< node 2 node projection
};

/// \enum BasisFunctions
/// \brief radial basis function weights for MMLS and SPLINE
enum BasisFunctions {
  WENDLAND2 = 0,  ///< Wendland 2nd order
  WENDLAND4,      ///< Wendland 4th order
  WENDLAND6,      ///< Wendland 6th order
  WU2,            ///< Wu 2nd order
  WU4,            ///< Wu 4th order
  WU6,            ///< Wu 6th order
  BUHMANN3        ///< Buhmann 3or order
};

/// \class Mapper
/// \brief the mapper interface for interface solution transfer
class Mapper {
  /// \brief initialize parameter list
  inline void init_parlist_() {
    FOR_DIR(i) {
      opts_[i].reset(new ::Teuchos::ParameterList());
      auto &list = *opts_[i];
      list.set("Map Type", "Point Cloud");
      auto &sub_list = list.sublist("Point Cloud", false);
      sub_list.set("Map Type", "Moving Least Square Reconstruction");
      sub_list.set("Spatial Dimension", dim_);
      sub_list.set("Basis Type", "Wendland");
      sub_list.set("Basis Order", 4);
      sub_list.set("Type of Search", "Radius");
      sub_list.set("RBF Radius", 0.0);
      sub_list.set("Num Neighbors", 0);
      sub_list.set("Matching Nodes", false);

      auto &sub_list_search = list.sublist("Search", false);
      sub_list_search.set("Track Missed Range Entities", true);
    }
  }

  /// \brief parse and formatting a parameter list
  /// \param[in] list parameter list
  inline static std::string parse_list_(::Teuchos::ParameterList &list) {
    using std::string;
    auto &sub_list = list.sublist("Point Cloud", true);
    const bool rs = sub_list.get<string>("Type of Search") == "Radius";
    std::ostringstream ss;
    ss << string(LEN2, ' ') << "map type: " << sub_list.get<string>("Map Type")
       << '\n'
       << string(LEN2, ' ')
       << "spatial dimension: " << sub_list.get<int>("Spatial Dimension")
       << '\n'
       << string(LEN2, ' ')
       << "basis type: " << sub_list.get<string>("Basis Type") << '\n'
       << string(LEN2, ' ')
       << "basis order: " << sub_list.get<int>("Basis Order") << '\n'
       << string(LEN2, ' ')
       << "type of search: " << sub_list.get<string>("Type of Search") << '\n'
       << string(LEN2, ' ') << (rs ? "radius: " : "knn: ")
       << (rs ? sub_list.get<double>("RBF Radius")
              : sub_list.get<int>("Num Neighbors"))
       << '\n';
    return ss.str();
  }

  /// \brief helper for set local search
  /// \tparam _Dir direction
  /// \tparam _V value type
  /// \param[in] type search type
  /// \param[in] value tag for value
  /// \param[in] v actual value
  template <bool _Dir, typename _V>
  inline void set_search(const std::string &type, const std::string &value,
                         const _V &v) noexcept {
    auto &list = opts_[_Dir]->sublist("Point Cloud", true);
    list.set("Type of Search", type);
    list.set(value, v);
  }

  /// \brief helper for get search info
  /// \tparam _Dir direction
  /// \tparam _V value type
  /// \param[in] type search type
  /// \param[in] value tag for value
  /// \param[in] dft default value
  template <bool _Dir, typename _V>
  inline _V get_search(const std::string &type, const std::string &value,
                       const _V &dft) const noexcept {
    const auto &list = opts_[_Dir]->sublist("Point Cloud", true);
    const std::string &tp = list.get<std::string>("Type of Search");
    if (tp != type) return dft;
    return list.get<_V>(value);
  }

 public:
  /// \brief constructor
  /// \param[in] B input blue mesh
  /// \param[in] G input green mesh
  /// \param[in] version passed from Python inteface
  /// \param[in] date passed from Python interface
  /// \param[in] profiling whether do simple profiling, i.e. wtime
  explicit Mapper(std::shared_ptr<IMeshDB> B, std::shared_ptr<IMeshDB> G,
                  const std::string &version = "", const std::string &date = "",
                  bool profiling = true)
      : B_(B),
        G_(G),
        dim_(3),
        ready_(false),
        profiling_(profiling),
        timer_(0.0) {
    // make sure the communicators are similar, both communicators should not
    // be null
    int flag;
    MPI_Comm_compare(B->comm(), G->comm(), &flag);
    throw_error_if(flag != MPI_IDENT,
                   "the green and blue comms must be identical");
    using std::string;
    // NOTE pass in time should have "%b %d %Y %H:%M:%S" format
    streamer_master(B_->rank())
        << '\n'
        << string(LEN1, '-') << "\n\n"
        << string(LEN2, ' ') << "parpydtk2 version: " << version << '\n'
        << string(LEN2, ' ') << "mapper created time: " << date << '\n'
        << string(LEN2, ' ') << "binary built time: " << __DATE__ << ' '
        << __TIME__ << '\n'
        << string(LEN2, ' ') << "total processes: " << B_->ranks() << "\n\n"
        << string(LEN1, '-') << '\n';

    // init par lists
    init_parlist_();
  }

  virtual ~Mapper() = default;

  /// \name mapper_py_interface
  ///@{

  /// \brief get the ranks
  inline int ranks() const noexcept { return B_->ranks(); }

  /// \brief get my rank
  inline int rank() const noexcept { return B_->rank(); }

  /// \brief get the communicator
  inline MPI_Comm comm() const noexcept { return B_->comm(); }

  /// \brief set dimension
  /// \param[in] dim geometry dimension
  inline void set_dimension(int dim) {
    B_->set_dimension(dim);
    G_->set_dimension(dim);
    dim_ = dim;
    FOR_DIR(i)
    opts_[i]->sublist("Point Cloud", true).set("Spatial Dimension", dim_);
  }

  /// \brief use moving least square, this is the default method
  /// \sa use_spline, use_n2n
  inline void use_mmls() noexcept {
    FOR_DIR(i)
    opts_[i]
        ->sublist("Point Cloud", true)
        .set("Map Type", "Moving Least Square Reconstruction");
  }

  /// \brief use spline interpolation method
  /// \sa use_mmls, use_n2n
  inline void use_spline() noexcept {
    FOR_DIR(i)
    opts_[i]
        ->sublist("Point Cloud", true)
        .set("Map Type", "Spline Interpolation");
  }

  /// \brief use node 2 node project
  /// \param[in] matching are the interfaces mathing?
  /// \sa use_mmls, use_spline
  inline void use_n2n(bool matching = false) noexcept {
    FOR_DIR(i) {
      auto &list = opts_[i]->sublist("Point Cloud", true);
      list.set("Map Type", "Node To Node");
      list.set("Matching Nodes", matching);
    }
  }

  /// \brief set basis function, default is Wendland 4th order
  /// \param[in] basis basis function and order
  /// \sa BasisFunctions
  inline void set_basis(int basis) {
    throw_error_if(basis < 0 || basis > BUHMANN3, "unknown method");
    if (basis != BUHMANN3) {
      const bool wld = basis < 3;
      const std::string bm = wld ? "Wendland" : "Wu";
      basis %= 3;
      const int order = 2 * (basis + 1);
      FOR_DIR(i) {
        auto &list = opts_[i]->sublist("Point Cloud", true);
        list.set("Basis Type", bm);
        list.set("Basis Order", order);
      }
      return;
    }
    FOR_DIR(i) {
      auto &list = opts_[i]->sublist("Point Cloud", true);
      list.set("Basis Type", "Buhmann");
      list.set("Basis Order", 3);
    }
  }

  /// \brief use knn for blue mesh
  /// \param[in] knn number of nearest neighbors
  inline void use_knn_b(int knn) {
    throw_error_if(knn <= 0, "invalid number of knn");
    show_warning_if(knn < 3, "potentially to small knn");
    set_search<true>("Nearest Neighbor", "Num Neighbors", knn);
  }

  /// \brief use knn for green mesh
  /// \param[in] knn number of nearest neighbors
  inline void use_knn_g(int knn) {
    throw_error_if(knn <= 0, "invalid number of knn");
    show_warning_if(knn < 3, "potentially to small knn");
    set_search<false>("Nearest Neighbor", "Num Neighbors", knn);
  }

  /// \brief use radius for blue
  /// \param[in] r physical domain radius support
  /// \sa use_knn
  inline void use_radius_b(double r) {
    throw_error_if(r <= 0.0, "invalid radius number");
    set_search<true>("Radius", "RBF Radius", r);
  }

  /// \brief use radius for green
  /// \param[in] r physical domain radius support
  /// \sa use_knn
  inline void use_radius_g(double r) {
    throw_error_if(r <= 0.0, "invalid radius number");
    set_search<false>("Radius", "RBF Radius", r);
  }

  /// \brief check method
  inline int check_method() const noexcept {
    // since we assign all parameters symmetrically, we just check one
    const std::string &method =
        opts_[0]->sublist("Point Cloud", true).get<std::string>("Map Type");
    if (method == "Moving Least Square Reconstruction") return MMLS;
    if (method == "Spline Interpolation") return SPLINE;
    return N2N;
  }

  /// \brief check basis
  inline int check_basis() const noexcept {
    const std::string &basis =
        opts_[0]->sublist("Point Cloud", true).get<std::string>("Basis Type");
    const int &order =
        opts_[0]->sublist("Point Cloud", true).get<int>("Basis Order");
    if (basis == "Wendland") switch (order) {
        case 2:
          return WENDLAND2;
        case 4:
          return WENDLAND4;
        case 6:
          return WENDLAND6;
      }

    if (basis == "Wu") switch (order) {
        case 2:
          return WU2;
        case 4:
          return WU4;
        case 6:
          return WU6;
      }
    return BUHMANN3;
  }

  /// \brief check blue knn
  /// \note if blue does not use knn, then negative value returned
  inline int knn_b() const noexcept {
    return get_search<true>("Nearest Neighbor", "Num Neighbors", -1);
  }

  /// \brief check green knn
  inline int knn_g() const noexcept {
    return get_search<false>("Nearest Neighbor", "Num Neighbors", -1);
  }

  /// \brief check blue radius
  /// \note if blue does not use radius, then -1.0 is returned
  inline double radius_b() const noexcept {
    return get_search<true>("Radius", "RBF Radius", -1.0);
  }

  /// \brief check green radius
  inline double radius_g() const noexcept {
    return get_search<false>("Radius", "RBF Radius", -1.0);
  }

  /// \brief get the dimension
  inline int dimension() const noexcept { return dim_; }

  /// \brief get blue mesh
  inline std::shared_ptr<IMeshDB> blue_mesh() const noexcept { return B_; }

  /// \brief get green mesh
  inline std::shared_ptr<IMeshDB> green_mesh() const noexcept { return G_; }

  /// \brief begin initialization
  inline void begin_initialization() {
    using std::string;
    throw_error_if(ready_, "this mapper has already been initialized");
    throw_error_if(!B_->ready() || !G_->ready(),
                   "at least one of the meshes is not ready");
    // parse parameter list
    string info[2];
    FOR_DIR(i)
    info[i] = parse_list_(*opts_[i]);
    streamer_master(B_->rank()) << '\n'
                                << string(LEN1, '-') << "\n\n"
                                << string(LEN2, ' ') << "blue ===> green:\n"
                                << info[1] << '\n'
                                << string(LEN2, ' ') << "green ===> blue:\n"
                                << info[0] << '\n'
                                << string(LEN1, '-') << '\n';
    ready_ = true;  // trigger flag here
    timer_ = 0.0;
  }

  /// \brief register coupling fields
  /// \param[in] bf blue meshdb field data
  /// \param[in] gf green meshdb field data
  /// \param[in] direct \a true for b->g, \a false for g->b
  inline void register_coupling_fields(const std::string &bf,
                                       const std::string &gf, bool direct) {
    throw_error_if(B_->field_dim(bf) != G_->field_dim(gf),
                   "field dimensions don\'t match");
    auto &b_dtk_fields = B_->dtk_fields();
    auto &g_dtk_fields = G_->dtk_fields();
    auto &src = direct ? b_dtk_fields[bf] : g_dtk_fields[gf];
    auto &tgt = direct ? g_dtk_fields[gf] : b_dtk_fields[bf];
    auto &optr = operators_[direct];
    auto build = optr.insert(
        std::make_pair(std::make_pair(bf, gf),
                       factory_.create(src.second->getMap(),
                                       tgt.second->getMap(), *opts_[direct])));
    if (!build.second) {
      show_warning(bf + "+" + gf + " already exists, ignoring request");
      return;
    }

    auto &op = build.first->second;
    auto &src_mngrs = direct ? B_->mangers() : G_->mangers();
    auto &tgt_mngrs = direct ? G_->mangers() : B_->mangers();
    auto &src_mngr = src_mngrs[src.first];
    auto &tgt_mngr = tgt_mngrs[tgt.first];

    // actually build the operator
    double t = MPI_Wtime();
    op->setup(src_mngr.functionSpace(), tgt_mngr.functionSpace());
    timer_ += MPI_Wtime() - t;
  }

  /// \brief check if a coupling data fields exists
  /// \param[in] bf blue meshdb field data
  /// \param[in] gf green meshdb field data
  /// \param[in] direct \a true for b->g, \a false for g->b
  inline bool has_coupling_fields(const std::string &bf, const std::string &gf,
                                  bool direct) {
    return operators_[direct].find(std::make_pair(bf, gf)) !=
           operators_[direct].end();
  }

  /// \brief end initialization
  inline void end_initialization() {
    throw_error_if(!ready_,
                   "the ready tag was not triggerred, did you forget "
                   "call begin_initialization?");
    double avg, min_, max_;
    if (profiling_ && ranks() > 1) {
      // gather all time information
      std::vector<double> ts(ranks());
      ts[0] = timer_;
      int ret = MPI_Gather(&timer_, 1, MPI_DOUBLE, ts.data(), 1, MPI_DOUBLE, 0,
                           comm());
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

      // compute timing data
      auto minmax = std::minmax_element(ts.cbegin(), ts.cend());
      min_ = *minmax.first;
      max_ = *minmax.second;
      avg = std::accumulate(ts.cbegin(), ts.cend(), 0.0) / ranks();
    } else if (profiling_)
      avg = min_ = max_ = timer_;

    // streaming messages
    using std::string;

    streamer_master(rank())
        << string(LEN1, '-') << "\n\n"
        << string(LEN2, ' ') << "total registered coupling fields: "
        << operators_[0].size() + operators_[1].size() << '\n'
        << string(LEN2, ' ') << "blue ===> green:\n";
    for (const auto &op : operators_[0])
      streamer_master(rank()) << string(LEN2 + 5, ' ') << op.first.first << ':'
                              << op.first.second << '\n';
    streamer_master(rank()) << string(LEN2, ' ') << "green ===> blue:\n";
    for (const auto &op : operators_[1])
      streamer_master(rank()) << string(LEN2 + 5, ' ') << op.first.second << ':'
                              << op.first.first << '\n';
    if (profiling_)
      streamer_master(rank())
          << string(LEN2, ' ') << "time used: " << std::scientific << "min "
          << min_ << ", max " << max_ << ", avg " << avg << '\n';
    streamer_master(rank()) << '\n' << string(LEN1, '-') << '\n';
  }

  /// \brief begin to transfer data
  inline void begin_transfer() noexcept { timer_ = 0.0; }

  /// \brief transfer data
  /// \param[in] bf blue meshdb field data
  /// \param[in] gf green meshdb field data
  /// \param[in] direct \a true for b->g, \a false for g->b
  inline void transfer_data(const std::string &bf, const std::string &gf,
                            bool direct) {
    const auto key = std::make_pair(bf, gf);
    auto op_iter = operators_[direct].find(key);
    throw_error_if(op_iter == operators_[direct].end(),
                   "no such operator exist");
    auto &b_dtk_fields = B_->dtk_fields();
    auto &g_dtk_fields = G_->dtk_fields();
    auto &src = direct ? b_dtk_fields[bf] : g_dtk_fields[gf];
    auto &tgt = direct ? g_dtk_fields[gf] : b_dtk_fields[bf];
    double t = MPI_Wtime();
    op_iter->second->apply(*src.second, *tgt.second);
    timer_ += MPI_Wtime() - t;
  }

  /// \brief end transfer
  inline void end_transfer() {
    double avg, min_, max_;
    if (profiling_ && ranks() > 1) {
      // gather all time information
      std::vector<double> ts(ranks());
      ts[0] = timer_;
      int ret = MPI_Gather(&timer_, 1, MPI_DOUBLE, ts.data(), 1, MPI_DOUBLE, 0,
                           comm());
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

      // compute timing data
      auto minmax = std::minmax_element(ts.cbegin(), ts.cend());
      min_ = *minmax.first;
      max_ = *minmax.second;
      avg = std::accumulate(ts.cbegin(), ts.cend(), 0.0) / ranks();
    } else if (profiling_)
      avg = min_ = max_ = timer_;

    using std::string;

    streamer_master(rank()) << string(LEN1, '-') << "\n\n"
                            << string(LEN2, ' ') << "transfer finished\n";
    if (profiling_)
      streamer_master(rank())
          << string(LEN2, ' ') << "time used: " << std::scientific << "min "
          << min_ << ", max " << max_ << ", avg " << avg << '\n';
    streamer_master(rank()) << '\n' << string(LEN1, '-') << '\n';
  }

  ///@}
 protected:
  /// \brief blue mesh
  std::shared_ptr<IMeshDB> B_;

  /// \brief green mesh
  std::shared_ptr<IMeshDB> G_;

  /// \brief dimension
  int dim_;

  /// \brief flag to indicate the mapper is ready for transfering
  bool ready_;

  /// \brief whether do simple profiling
  bool profiling_;

  /// \brief a simple timer buffer
  double timer_;

  /// \brief parameter list
  std::unique_ptr<Teuchos::ParameterList> opts_[2];

  /// \brief transfer operators
  std::map<std::pair<std::string, std::string>,
           Teuchos::RCP<DataTransferKit::MapOperator>>
      operators_[2];

  /// \brief map factory
  static ::DataTransferKit::MapOperatorFactory factory_;
};

// define the factory
::DataTransferKit::MapOperatorFactory Mapper::factory_;

/** @}*/

}  // namespace parpydtk2

#undef FOR_DIR
