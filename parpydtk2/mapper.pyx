#!python
#cython: language_level=3, boundscheck=False, embedsignature=True, wraparound=False

"""DTK2 mapper interface with MOAB"""

# cimports
cimport numpy as cnp
from libcpp cimport bool
from libcpp.vector cimport vector as std_vector
from libcpp.string cimport string as std_string
from mpi4py cimport MPI as c_MPI
from mpi4py cimport libmpi as mpi

ctypedef c_MPI.Comm comm_t
ctypedef mpi.MPI_Comm c_comm_t


cdef extern from 'mpi-compat.hpp':
    pass

# normal imports
import numpy as np
import mpi4py
import datetime
from ._version import __version__


__version__ = __version__
__author__ = 'Qiao Chen'
__copyright__ = 'Copyright 2018, Qiao Chen'


cdef extern from 'src/dtk2.hpp' namespace 'parpydtk2' nogil:


    cdef cppclass _IMeshDB 'parpydtk2::IMeshDB':
        _IMeshDB(c_comm_t comm) except +
        void begin_create() except +
        void create_vset() except +
        void create_vertices(int nv, const double *coords, unsigned set_id) except +
        void assign_gids(int nv, const int *gids, unsigned set_id) except +
        void finish_create(bool trivial_gid) except +
        int size(unsigned set_id) except +
        int sets()
        void get_bbox(double *v, unsigned set_id) except +
        void get_gbbox(double *v, unsigned set_id) except +
        void create_field(const std_string &field_name, unsigned set_id, int dim) except +
        bool has_field(const std_string &field_name)
        int field_dim(const std_string &field_name) except +
        int field_set_id(const std_string &field_name) except +
        void assign_field(const std_string &field_name, const double *values,
            unsigned set_id) except +
        void extract_field(const std_string &field_name, double *values,
            unsigned set_id) except +


    cdef cppclass _Mapper 'parpydtk2::Mapper':
        _Mapper(c_comm_t comm, const std_string &version,
            const std_string &date, bool profiling) except +
        int ranks()
        int rank()
        c_comm_t comm()
        void set_dimension(int dim) except +
        void use_mmls()
        void use_spline()
        void use_n2n(bool matching)
        void set_basis(int basis) except +
        void use_knn_b(int knn) except +
        void use_knn_g(int knn) except +
        void use_radius_b(double r) except +
        void use_radius_g(double r) except +
        int check_method()
        int check_basis()
        int knn_b()
        int knn_g()
        double radius_b()
        double radius_g()
        int dimension()
        _IMeshDB &blue_mesh()
        _IMeshDB &green_mesh()
        void begin_initialization() except +
        void register_coupling_fields(const std_string &bf,
            const std_string &gf, bool direct) except +
        bool has_coupling_fields(const std_string &bf,
            const std_string &gf, bool direct)
        void end_initialization() except +
        void begin_transfer()
        void transfer_data(const std_string &bf,
            const std_string &gf, bool direct) except +
        void end_transfer()


cdef class IMeshDB:
    cdef _IMeshDB *mdb
    cdef comm_t comm

    # dummy constructor for doc
    def __init__(self, comm=None):
        pass

    def __cinit__(self, comm_t comm=None):
        cdef c_comm_t comm_ = mpi.MPI_COMM_WORLD if comm is None else comm.ob_mpi
        self.comm = c_MPI.Comm()
        self.comm.ob_mpi = comm_
        self.mdb = new _IMeshDB(comm_)

    def __dealloc__(self):
        del self.mdb

    @property
    def comm(self):
        """MPI.Comm: communicator"""
        return self.comm

    def begin_create(self):
        """Begin to create/manupilate the mesh

        This function must be called in order to let the mesh databse be
        aware that you will create meshes.

        See Also
        --------
        :func:`finish_create` : finish creating mesh
        """
        self.mdb.begin_create()

    def create_vset(self):
        """Create a new vertex set

        .. note:: the default set is root set

        See Also
        --------
        :func:`create_vertices` : create vertices
        """
        self.mdb.create_vset()

    def create_vertices(self, double[:, ::1] coords, unsigned set_id=0):
        """Create a set of coordinates

        Parameters
        ----------
        coords : np.ndarray
            nx3 coordinates in double precision
        set_id : int
            set id, default is the root set

        See Also
        --------
        :func:`create_vset` : create vertex set
        :func:`assign_gids` : assign global IDs
        """
        assert coords.shape[1] == 3
        self.mdb.create_vertices(
            <int> coords.shape[0],
            <const double *> &coords[0, 0],
            set_id
        )

    def assign_gids(self, int[::1] gids, unsigned set_id=0):
        """Assign global IDs

        .. todo:: enrich doc

        Parameters
        ----------
        gids : np.ndarray
            global IDs
        set_id : int
            set id, default is the root set

        See Also
        --------
        :func:`create_vertices` : create vertices
        """
        self.mdb.assign_gids(
            <int> gids.size,
            <const int *> &gids[0],
            set_id
        )

    def finish_create(self, trivial_gid=True):
        """finish mesh creation

        .. todo:: enrich doc

        Parameters
        ----------
        trivial_gid : bool
            `True` if we use MOAB trivial global ID computation
        """
        cdef bool tg = <bool> 1 if trivial_gid else <bool> 0
        self.mdb.finish_create(tg)

    def size(self, unsigned set_id=0):
        """Get the size of a set

        Parameters
        ----------
        set_id : int
            index os set id

        Returns
        -------
        int
            number of vertices in set `set_id`

        See Also
        --------
        :func:`sets` : get the number of sets
        """
        return self.mdb.size(set_id)

    @property
    def sets(self):
        """int: Get the number of sets"""
        return self.mdb.sets()

    def bbox(self, unsigned set_id=0):
        """Get the bounding box [[min1,min2,min3],[max1,max2,max3]]

        Parameters
        ----------
        set_id : int
            set index

        Returns
        -------
        np.ndarray
            2x3 bounding box
        """
        cdef cnp.ndarray[double, ndim=2, mode='c'] box = \
            np.empty(shape=(2, 3), dtype='double')
        self.mdb.get_bbox(<double *> box.data, set_id)
        return box

    def gbbox(self, unsigned set_id=0):
        """Get the global bounding box [[min1,min2,min3],[max1,max2,max3]]

        Parameters
        ----------
        set_id : int
            set index

        Returns
        -------
        np.ndarray
            2x3 bounding box
        """
        cdef cnp.ndarray[double, ndim=2, mode='c'] box = \
            np.empty(shape=(2, 3), dtype='double')
        self.mdb.get_gbbox(<double *> box.data, set_id)
        return box

    def create_field(self, str field_name, unsigned set_id=0, int dim=1):
        """Create a data field for solution transfer

        Parameters
        ----------
        field_name : str
            name of the field
        set_id : int
            index os set id
        dim : int
            dimension of the field, i.e. scalar, vector, tensor
        """
        cdef std_string fn = <std_string> field_name.encode('UTF-8')
        self.mdb.create_field(fn, set_id, dim)

    def has_field(self, str field_name):
        """Check if a field exists

        Parameters
        ----------
        field_name : str
            name of the field

        Returns
        -------
        bool
            `True` if this meshdb has `field_name`
        """
        cdef std_string fn = <std_string> field_name.encode('UTF-8')
        return self.mdb.has_field(fn)

    def field_dim(self, str field_name):
        """Check the field data dimension

        Parameters
        ----------
        field_name : str
            name of the field

        Returns
        -------
        int
            data field dimension of `field_name`
        """
        cdef std_string fn = <std_string> field_name.encode('UTF-8')
        return self.mdb.field_dim(fn)

    def field_set_id(self, str field_name):
        """Check the field set id

        Parameters
        ----------
        field_name : str
            name of the field

        Returns
        -------
        int
            data field set id of `field_name`
        """
        cdef std_string fn = <std_string> field_name.encode('UTF-8')
        return self.mdb.field_set_id(fn)

    def assign_field(self, str field_name, cnp.ndarray values not None,
        unsigned set_id=0):
        """Assign values to a field

        .. note:: `values` size must be at least size*dim

        Parameters
        ----------
        field_name : str
            name of the field
        values : np.ndarray
            input source values
        set_id : int
            set index, default is root set

        See Also
        --------
        :func:`extract_field` : extract value from a field
        :func:`size` : check the size of a mesh set
        :func:`dim` : get the dimension
        """
        cdef:
            std_string fn = <std_string> field_name.encode('UTF-8')
            int n = self.mdb.size(set_id) * self.mdb.field_dim(fn)
        assert values.size >= n
        self.mdb.assign_field(fn, <const double *> values.data, set_id)

    def extract_field(self, str field_name, unsigned set_id=0,
        double[::1] buffer=None, reshape=False):
        """Extact the values from a field

        .. warning:: if `buffer` is passed in, it must be 1D

        Parameters
        ----------
        field_name : str
            name of the field
        set_id : int
            set index, default is root set
        buffer : np.ndarray
            1D buffer
        reshape : bool
            `True` if we reshape the output, only for vectors/tensors

        Returns
        -------
        np.ndarray
            field data values
        """
        cdef:
            std_string fn = <std_string> field_name.encode('UTF-8')
            int n = self.mdb.size(set_id)
            int dim = self.mdb.field_dim(fn)
            cnp.ndarray[double, ndim=1] values
        if buffer is not None:
            assert len(buffer) >= n * dim
            values = np.asarray(buffer)
        else:
            values = np.empty(n * dim, dtype='double')
        self.mdb.extract_field(fn, <double *> values.data, set_id)
        if reshape and dim > 1:
            return values.reshape((n, dim))
        return values


cdef class Mapper:
    cdef _Mapper *mp

    # dummpy contructor for doc
    def __init__(self, comm=None, profiling=True):
        pass

    def __cinit__(self, comm_t comm=None, profiling=True):
        cdef:
            std_string version = __version__.encode('UTF-8')
            std_string date = \
                datetime.datetime.now().strftime('%b %d %Y %H:%M:%S').encode('UTF-8')
            bool prof = <bool> 1 if profiling else <bool> 0
            c_comm_t comm_
        if comm is None:
            comm_ = mpi.MPI_COMM_WORLD
        else:
            comm_ = comm.ob_mpi
        self.mp = new _Mapper(comm_, version, date, prof)

    def __dealloc__(self):
        del self.mp

    @property
    def ranks(self):
        """int: Check the total process number

        See Also
        --------
        :attr:`rank` : "my" rank
        """
        return self.mp.ranks()

    @property
    def rank(self):
        """int: Check "my" rank

        See Also
        --------
        :attr:`ranks` : get the total communicator size
        :attr:`comm` : MPI communicator
        """
        return self.mp.rank()

    @property
    def comm(self):
        """MPI.Comm: Get the communicator

        See Also
        --------
        :attr:`ranks` : get the total communicator size
        :attr:`rank` : "my" rank
        """
        cdef c_comm_t comm_ = self.mp.comm()
        cdef comm_t comm = c_MPI.Comm()
        comm.ob_mpi = comm_
        return comm

    @property
    def dimension(self):
        """int: Get the problem dimension

        .. note:: this is the spacial dimension
        """
        return self.mp.dimension()

    @dimension.setter
    def dimension(self, int dim):
        self.mp.set_dimension(dim)

    @property
    def blue_mesh(self):
        """:class:`IMeshDB`: Get the blue side mesh database

        .. note:: this call is where you can access the meshdb

        See Also
        --------
        :attr:`green_mesh` : get the green mesh
        """
        cdef IMeshDB mdb = IMeshDB()
        mdb.mdb = <_IMeshDB *> &self.mp.blue_mesh()
        return mdb

    @property
    def green_mesh(self):
        """:class:`IMeshDB`: Get the blue side mesh database

        .. note:: this call is where you can access the meshdb

        See Also
        --------
        :attr:`blue_mesh` : get the blue mesh
        """
        cdef IMeshDB mdb = IMeshDB()
        mdb.mdb = <_IMeshDB *> &self.mp.green_mesh()
        return mdb

    @property
    def method(self):
        """int: Get the method tag

        See Also
        --------
        :attr:`basis` : the basis function and order attribute
        """
        return self.mp.check_method()

    @method.setter
    def method(self, int m):
        if m == 0:
            # mls
            self.mp.use_mmls()
        elif m == 1:
            self.mp.use_spline()
        elif m == 2:
            self.mp.use_n2n(<bool> 0)
        else:
            raise AttributeError('unknown method')

    def set_matching_flag_n2n(self, bool matching):
        """Set the matching flag for N2N

        .. note:: this function will not throw even if you dont use n2n

        Parameters
        ----------
        matching : bool
            `True` if the interfaces are matching
        """
        cdef bool flag = <bool> 1 if matching else <bool> 0
        self.mp.use_n2n(flag)

    @property
    def basis(self):
        """int: Get the basis function flag

        See Also
        --------
        :attr:`method` : get the method tag
        """
        return self.mp.check_basis()

    @basis.setter
    def basis(self, int bf):
        self.mp.set_basis(bf)

    @property
    def knn_b(self):
        """int: KNN of blue mesh

        .. note:: if blue does not use KNN, then -1 returned

        See Also
        --------
        :attr:`blue_mesh` : blue mesh database
        :attr:`knn_g` : green knn
        """
        return self.mp.knn_b()

    @knn_b.setter
    def knn_b(self, int knn):
        self.mp.use_knn_b(knn)

    @property
    def knn_g(self):
        """int: KNN of green mesh

        .. note:: if green does not use KNN, then -1 returned

        See Also
        --------
        :attr:`green_mesh` : green mesh database
        :attr:`knn_b` : blue knn
        """
        return self.mp.knn_g()

    @knn_g.setter
    def knn_g(self, int knn):
        self.mp.use_knn_g(knn)

    @property
    def radius_b(self):
        """float: physical domain radius support for blue mesh

        .. note:: if blue does not use RBF-search, then -1.0 returned

        See Also
        --------
        :attr:`blue_mesh` : blue mesh database
        :attr:`radius_g` : green radius
        """
        return self.mp.radius_b()

    @radius_b.setter
    def radius_b(self, double r):
        self.mp.use_radius_b(r)

    @property
    def radius_g(self):
        """float: physical domain radius support for green mesh

        .. note:: if green does not use RBF-search, then -1.0 returned

        See Also
        --------
        :attr:`green_mesh` : green mesh database
        :attr:`radius_b` : blue radius
        """
        return self.mp.radius_g()

    @radius_g.setter
    def radius_g(self, double r):
        self.mp.use_radius_g(r)

    def begin_initialization(self):
        """Initialization starter

        This is a must-call function in order to indicate mapper that you are
        about to initialize/register coupling fields

        See Also
        --------
        :func:`register_coupling_fields` : register coupled fields
        :func:`end_initialization` : finish initialization
        """
        self.mp.begin_initialization()

    def register_coupling_fields(self, str bf, str gf, bool direct):
        """register a coupled fields

        .. note:: we use boolean to indicate direction

        Parameters
        ----------
        bf : str
            blue mesh field name
        gf : str
            green mesh field name
        direct : bool
            `True` for blue->green, `False` for the opposite

        See Also
        --------
        :func:`transfer_data` : transfer a coupled data fields
        """
        cdef:
            std_string bf_ = <std_string> bf.encode('UTF-8')
            std_string gf_ = <std_string> gf.encode('UTF-8')
            bool dr = <bool> 1 if direct else <bool> 0
        self.mp.register_coupling_fields(bf_, gf_, dr)

    def has_coupling_fields(self, str bf, str gf, bool direct):
        """Check if a coupled fields exists

        Returns
        -------
        bool
            `True` if (bf,gf) exists in `direct`ion
        """
        cdef:
            std_string bf_ = <std_string> bf.encode('UTF-8')
            std_string gf_ = <std_string> gf.encode('UTF-8')
            bool dr = <bool> 1 if direct else <bool> 0
        return self.mp.has_coupling_fields(bf_, gf_, dr)

    def end_initialization(self):
        """Initialization closer

        This is a must-call function in order to tell the mapper we are ready

        See Also
        --------
        :func:`begin_initialization` : initialization starter
        """
        self.mp.end_initialization()

    def begin_transfer(self):
        """Transfer starter

        This is a must-call function to inidate the beginning of a transferring
        block

        See Also
        --------
        :func:`end_transfer` : transfer closer
        :func:`tranfer_data` : transfer a coupled data fields
        """
        self.mp.begin_transfer()

    def transfer_data(self, str bf, str gf, bool direct):
        """Transfer (bf, gf) in `direct`ion

        Parameters
        ----------
        bf : str
            blue mesh field name
        gf : str
            green mesh field name
        direct : bool
            `True` for blue->green, `False` for the opposite

        See Also
        --------
        :func:`register_coupling_fields` : register coupled fields
        """
        cdef:
            std_string bf_ = <std_string> bf.encode('UTF-8')
            std_string gf_ = <std_string> gf.encode('UTF-8')
            bool dr = <bool> 1 if direct else <bool> 0
        self.mp.transfer_data(bf_, gf_, dr)

    def end_transfer(self):
        """Transfer closer

        This is a must-call function to indicate we have finished a sequence of
        transferring requests

        See Also
        --------
        :func:`begin_transfer` : transfer starter
        """
        self.mp.end_transfer()
