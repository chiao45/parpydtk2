#ifndef _PYDTK2_DTKWRAP_HPP
#define _PYDTK2_DTKWRAP_HPP

#include <cassert>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <moab/Core.hpp>
#include <moab/Interface.hpp>
#include <moab/ParallelComm.hpp>

#include <DTK_FunctionSpace.hpp>
#include <DTK_MapOperatorFactory.hpp>
#include <DTK_MoabManager.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_MultiVector.hpp>

// adjust a moab macro
#undef MB_CHK_ERR
#define MB_CHK_ERR(err_code)                                                   \
  do {                                                                         \
    if (moab::MB_SUCCESS != err_code)                                          \
      std::exit(moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__,    \
                              err_code, "", moab::MB_ERROR_TYPE_EXISTING));    \
  } while (false)

namespace pydtk2 {

// This is the header for warpping DTK with MOAB serial version
// we employ the Moving Least Square with radial basis functions
// method for effieciency

// since DTK only works in "parallel", check MPI is required. However,
// a mapper should not be responsible for MPI in general.
static int mpi_flag = 0, mpi_init = 0;

#define CHECK_INIT_MPI                                                         \
  mpi_init = MPI_Initialized(&mpi_flag);                                       \
  if (mpi_init != MPI_SUCCESS || !mpi_flag) {                                  \
    int argc = 0;                                                              \
    char **argv = nullptr;                                                     \
    MPI_Init(&argc, &argv);                                                    \
  }

class Dtk2MoabManager {
public:
  // constructors
  Dtk2MoabManager() = delete;

  // real constructors
  Dtk2MoabManager(::moab::Interface *meshdb, unsigned long mesh_set = 0)
      : _meshdb(meshdb) {
    ::moab::ErrorCode err;

    // configure the mesh set
    if (mesh_set)
      _mesh_set = mesh_set;
    else
      _mesh_set = _meshdb->get_root_set();

    // retrieve vertices
    err = _meshdb->get_entities_by_type(_mesh_set, ::moab::MBVERTEX, _nodes);
    MB_CHK_ERR(err);

    // create a meshset only contains the nodes
    err = _meshdb->create_meshset(0, _node_set);
    MB_CHK_ERR(err);
    err = _meshdb->add_entities(_node_set, _nodes.data(), _nodes.size());
    MB_CHK_ERR(err);

    // check mpi
    CHECK_INIT_MPI;

    // create the dummy parallel moab manager
    _pmeshdb.reset(
        new (std::nothrow)::moab::ParallelComm(_meshdb, MPI_COMM_WORLD));
    Teuchos::RCP<::moab::ParallelComm> mdb =
        Teuchos::rcp(::moab::ParallelComm::get_pcomm(_meshdb, 0), false);
    _dtk_moab_manager.reset(new (std::nothrow)
                                DataTransferKit::MoabManager(mdb, _mesh_set));

    // reserve spaces for vector fuctions, and operators
    _v_pool.reserve(20);
  }

  // register field associated with tag in moab by passing in the string value
  // of the tags
  void register_tag(const std::string &tag) {
    ::moab::ErrorCode err;
    ::moab::Tag mtag;

    // get the moab base from par moab, then retrieve the tag
    err = _pmeshdb->get_moab()->tag_get_handle(tag.c_str(), mtag);
    MB_CHK_ERR(err);

    // inset to the set of tags strings
    if (_tags.find(tag) != _tags.end())
      return;
    _tags[tag] = _v_pool.size();

    // create the dtk vector function space accordingly
    _v_pool.push_back(
        _dtk_moab_manager->createFieldMultiVector(_node_set, mtag));
    _tag_sizes.push_back(mtag->get_size() /
                         mtag->size_from_data_type(mtag->get_data_type()));
  }

  // utility checking tags
  bool has_tag(const std::string &tag) const noexcept {
    return _tags.find(tag) != _tags.end();
  }

private:
  // MOAB base
  ::moab::Interface *_meshdb;
  // parallel managers, dummy for now
  std::unique_ptr<::moab::ParallelComm> _pmeshdb;
  // tags
  std::unordered_map<std::string, unsigned> _tags;
  // mesh set, and node set
  ::moab::EntityHandle _mesh_set, _node_set;
  // dtk2moab manager, and dtk vector spaces
  std::unique_ptr<::DataTransferKit::MoabManager> _dtk_moab_manager;
  std::vector<::Teuchos::RCP<::DataTransferKit::FieldMultiVector>> _v_pool;
  std::vector<int> _tag_sizes;

  // store nodes entity handles
  std::vector<::moab::EntityHandle> _nodes;

  friend class Dtk2Mapper;
};

// The actual mapper class

static DataTransferKit::MapOperatorFactory factory;

class Dtk2Mapper {
public:
  // constructors
  Dtk2Mapper() = delete;

  // real constructor
  // for moving least square, there are 3 basis functions
  //  wenland := 0; wu := 1; and buhmann := 2
  // the acceptable orders for each one are, {0, 2, 4, 6}, {0, 2, 4}, {3}, resp
  Dtk2Mapper(Dtk2MoabManager *source, Dtk2MoabManager *target, double r,
             int basis = 0, int order = 4, int dim = 3,
             bool track_missing = true)
      : _source(source), _target(target),
        _options(new (std::nothrow)::Teuchos::ParameterList()) {
    ::moab::ErrorCode err;
    err = _source->_pmeshdb->get_moab()->set_dimension(dim);
    MB_CHK_ERR(err);
    err = _target->_pmeshdb->get_moab()->set_dimension(dim);
    MB_CHK_ERR(err);

    // Basis and order debuggings are done in python scope!
    std::string basis_string = Dtk2Mapper::_get_basis_from_id(basis);

    // configurations
    // what we use is moving least square with radial basis functions.
    _options->set("Map Type", "Point Cloud");
    auto &sub_list = _options->sublist("Point Cloud", false);
    sub_list.set("Map Type", "Moving Least Square Reconstruction");
    sub_list.set("Spatial Dimension", dim);
    sub_list.set("Basis Type", basis_string);
    sub_list.set("Basis Order", order);
    sub_list.set("Type of Search", "Radius");
    sub_list.set("RBF Radius", r);

    auto &sub_list_search = _options->sublist("Search", false);
    sub_list_search.set("Track Missed Range Entities", track_missing);

    _op_pool.reserve(20);
  }

  // set methods
  void set_basis(int basis) {
    std::string basis_string = Dtk2Mapper::_get_basis_from_id(basis);
    auto &sub_list = _options->sublist("Point Cloud", true);
    sub_list.set("Basis Type", basis_string);
  }

  void set_spacial_dimension(int dim) {
    auto &sub_list = _options->sublist("Point Cloud", true);
    sub_list.set("Spatial Dimension", dim);
  }

  void set_order(int order) {
    auto &sub_list = _options->sublist("Point Cloud", true);
    sub_list.set("Basis Order", order);
  }

  void set_rbf_radius(double r) {
    auto &sub_list = _options->sublist("Point Cloud", false);
    sub_list.set("RBF Radius", r);
  }

  void set_track_missing_flag(bool flag) {
    auto &sub_list_search = _options->sublist("Search", true);
    sub_list_search.set("Track Missed Range Entities", flag);
  }

  // register coupled tags
  void register_coupled_tags(const std::string &stag, const std::string &ttag) {
    assert(_source->has_tag(stag) && "Unknown source tag.");
    assert(_target->has_tag(ttag) && "Unknown target tag.");

    // check if already registered
    if (has_coupled_tags(stag, ttag))
      return;
    _coupled_tags[std::make_pair(stag, ttag)] = _op_pool.size();
    const unsigned skey = _get_key(stag, true);
    const unsigned tkey = _get_key(ttag, false);

    // set up operator
    _op_pool.push_back(factory.create(_source->_v_pool[skey]->getMap(),
                                      _target->_v_pool[tkey]->getMap(),
                                      *_options));
    _op_pool.back()->setup(_source->_dtk_moab_manager->functionSpace(),
                           _target->_dtk_moab_manager->functionSpace());
  }

  // transfer the values
  void apply(const std::string &stag, const std::string &ttag) {
    assert(has_coupled_tags(stag, ttag) && "Unknown coupled tags.");

    // the followings may be inefficient
    const unsigned key = _coupled_tags[std::make_pair(stag, ttag)];
    const unsigned skey = _get_key(stag, true);
    const unsigned tkey = _get_key(ttag, false);

    _op_pool[key]->apply(*_source->_v_pool[skey], *_target->_v_pool[tkey]);
  }

  bool has_coupled_tags(const std::string &stag, const std::string &ttag) {
    return _coupled_tags.find(std::make_pair(stag, ttag)) !=
           _coupled_tags.end();
  }

private:
  Dtk2MoabManager *_source, *_target;
  // configurations
  std::unique_ptr<::Teuchos::ParameterList> _options;
  // transfer oprators
  typedef std::pair<std::string, std::string> coupled_key_t;
  std::map<coupled_key_t, int> _coupled_tags;
  std::vector<::Teuchos::RCP<DataTransferKit::MapOperator>> _op_pool;

  static std::string _get_basis_from_id(int basis);

  unsigned _get_key(const std::string &tag, bool source) const noexcept {
    if (source)
      return _source->_tags[tag];
    else
      return _target->_tags[tag];
  }
};

std::string Dtk2Mapper::_get_basis_from_id(int basis) {
  return basis == 0 ? "Wendland"
                    : basis == 1 ? "Wu" : basis == 2 ? "Buhmann" : "";
}

} // namespace pydtk2

#endif
