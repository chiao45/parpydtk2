#!python
#cython: language_level=3, boundscheck=False, embedsignature=True, wraparound=False

"""DTK2 mapper interface with MOAB"""

# cimports
cimport numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector as std_vector
from libcpp.string cimport string as std_string
from mpi4py cimport MPI

# normal imports
import sys
import numpy as np
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI
import datetime
from ._version import __version__


cdef extern from 'src/dtk2.hpp' namespace 'parpydtk2':


    cdef cppclass _IMeshDB 'IMeshDB':
        void begin_create() except +
        void create_vset() except +
        void create_vertices(int nv, const double *coords, unsigned set_id) except +
        void assign_gids(int nv, const int *gids, unsigned set_id) except +
        void finish_create(bool trivial_gid) except +
        int size(unsigned set_id) except +
        int sets()
        void create_field(const std_string &field_name, unsigned set_id set, int dim) except +
        bool has_field(const std_string &field_name)
        int field_dim(const std_string &field_name) except +
        int field_set_id(const std_string &field_name) except +
        void assign_field(const std_string &field_name, const double *values,
            unsigned set_id) except +
        void extract_field(const std_string &field_name, double *values,
            unsigned set_id) except +


    cdef enum Methods:
        MMLS = 0
        SPLINE
        N2N


    cdef enum BasisFunctions:
        WENDLAND2 = 0
        WENDLAND4
        WENDLAND6
        WU2
        WU4
        WU6
        BUHMANN3


    cdef cppclass _Mapper 'Mapper':
        Mapper(MPI.MPI_Comm comm, const std_string &version,
            const std_string &date, bool profiling) except +
        int ranks()
        int rank()
        MPI.MPI_Comm comm()
        void set_dimension(int dim) except +
        _IMeshDB &blue_mesh()
        _IMeshDB &green_mesh()
        void begin_initialization() except +
        void register_coupling_fields(const std_string &bf,
            const std_string &gf, bool direct) except +
        void end_initialization() except +
        void begin_transfer()
        void tranfer_data(const std_string &bf,
            const std_string &gf, bool direct) except +
        void end_transfer()


MMLS = Methods.MMLS
SPLINE = Methods.SPLINE
N2N = Methods.N2N


WENDLAND2 = BasisFunctions.WENDLAND2
WENDLAND4 = BasisFunctions.WENDLAND4
WENDLAND6 = BasisFunctions.WENDLAND6
WU2 = BasisFunctions.WU2
WU4 = BasisFunctions.WU4
WU6 = BasisFunctions.WU6
BUHMANN3 = BasisFunctions.BUHMANN3


cdef class IMeshDB:
    cdef _IMeshDB *mdb

    # dummy constructor for doc
    def __init__(self):
        pass

    def __cinit__(self):
        self.mdb = NULL

    def __dealloc__(self):
        pass

    def begin_create(self):
        """Begin to create/manupilate the mesh

        This function must be called in order to let the mesh databse be
        aware that you will create meshes.

        See Also
        --------
        :func:`finish_create` : finish creating mesh
        """
        try:
            self.mdb.begin_create()
        except Exception as e:
            if MPI.Is_initialized():
                print('fail to create mesh: %s' % e, file=sys.stderr)
                MPI.COMM_WORLD.Abort()
            else:
                raise

    def create_vset(self):
        """Create a new vertex set

        .. note:: the default set is root set

        See Also
        --------
        :func:`create_vertices` : create vertices
        """
        try:
            self.mdb.create_vset()
        except Exception as e:
            if MPI.Is_initialized():
                print('fail to create vset: %s' % e, file=sys.stderr)
                MPI.COMM_WORLD.Abort()
            else:
                raise

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
        try:
            self.mdb.create_vertices(
                <int> coords.shape[0],
                <const double *> &coords[0],
                set_id
            )
        except Exception as e:
            if MPI.Is_initialized():
                print('fail to create vertices: %s' % e, file=sys.stderr)
                MPI.COMM_WORLD.Abort()
            else:
                raise

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
        try:
            self.mdb.create_gids(
                <int> gids.size,
                <const int *> &coords[0],
                set_id
            )
        except Exception as e:
            if MPI.Is_initialized():
                print('fail to create gids: %s' % e, file=sys.stderr)
                MPI.COMM_WORLD.Abort()
            else:
                raise

    def finish_create(self, trivial_gid=True):
        """finish mesh creation

        .. todo:: enrich doc

        Parameters
        ----------
        trivial_gid : bool
            `True` if we use MOAB trivial global ID computation
        """
        cdef bool tg = <bool> 1 if trivial_gid else <bool> 0
        try:
            self.mdb.finished(tg)
        except Exception as e:
            if MPI.Is_initialized():
                print('fail to finish create: %s' % e, file=sys.stderr)
                MPI.COMM_WORLD.Abort()
            else:
                raise

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
        try:
            return self.mdb.size(set_id)
        except Exception as e:
            if MPI.Is_initialized():
                print('fail to get size: %s' % e, file=sys.stderr)
                MPI.COMM_WORLD.Abort()
            else:
                raise

    def sets(self):
        """Get the nunber of sets

        Returns
        -------
        int
            number of sets
        """
        return self.mdb.sets()

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
        try:
            self.mdb.create_field(fn, set_id, dim)
        except Exception as e:
            if MPI.Is_initialized():
                print('fail to create field: %s' % e, file=sys.stderr)
                MPI.COMM_WORLD.Abort()
            else:
                raise

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


cdef class Mapper:
    cdef _Mapper *mp

    # dummpy contructor for doc
    def __init__(self, comm=MPI.COMM_WORLD, profiling=True):
        pass

    def __cinit__(self, MPI.MPI_Comm comm=MPI.COMM_WORLD, profiling=True):
        cdef:
            std_string version = __version__.decode('UTF-8')
            std_string date = \
                datetime.datetime.now().strftime('%b %d %Y %H:%M:%S').decode('UTF-8')
            bool prof = <bool> 1 if profiling else <bool> 0
        try:
            self.mp = new _Mapper(comm, version, date, prof)
        except Exception as e:
            if MPI.Is_initialized():
                print('fail to constructor mapper: %s' % e, file=sys.stderr)
                MPI.COMM_WORLD.Abort()
            else:
                raise

    def __dealloc__(self):
        del self.mp

    def ranks(self):
        """Check the ranks

        Returns
        -------
        int
            total process number
        """
        return self.mp.ranks()

    def rank(self):
        """Check "my" rank

        Returns
        -------
        int
            my rank

        See Also
        --------
        :func:`ranks` : get the total communicator size
        """
        return self.mp.rank()

    def comm(self):
        """Get the communicator

        Returns
        -------
        MPI.Comm
            MPI communicator
        """
        cdef MPI.MPI_Comm comm = self.mp.comm()
        return comm

    def set_dimension(self, int dim):
        """Set the spacial dimension

        .. note:: the default value is 3
        """
        try:
            self.mp.set_dimension(dim)
        except Exception as e:
            if MPI.Is_initialized():
                print('fail to set dimension: %s' % e, file=sys.stderr)
                MPI.COMM_WORLD.Abort()
            else:
                raise

    @property
    def blue_mesh(self):
        """Get the blue side mesh

        .. note:: this call is where you can access the meshdb

        Returns
        -------
        :class:`IMeshDB`
            mesh database of blue side

        See Also
        --------
        :attr:`green_mesh` : get the green mesh
        """
        cdef IMeshDB mdb
        mdb.mdb = <_IMeshDB *> &self.mp.blue_mesh()
        return mdb

    @property
    def green_mesh(self):
        """Get the blue side mesh

        .. note:: this call is where you can access the meshdb

        Returns
        -------
        :class:`IMeshDB`
            mesh database of blue side

        See Also
        --------
        :attr:`blue_mesh` : get the blue mesh
        """
        cdef IMeshDB mdb
        mdb.mdb = <_IMeshDB *> &self.mp.green_mesh()
        return mdb

    def begin_initialization(self):
        pass
