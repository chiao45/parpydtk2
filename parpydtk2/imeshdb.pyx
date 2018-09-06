#!python
#cython: language_level=3
#cython: boundscheck=False
#cython: embedsignature=True
#cython: wraparound=False

"""Interface mesh database

This module defines the representation of an interface mesh databse patch. It
is built on top of MOAB parallel mesh database. Since we are only interested in
the meshless methods, only point clouds are needed in :class:`IMeshDB`.

.. moduleauthor:: Qiao Chen <benechiao@gmail.com>
"""

cimport numpy as cnp
from cython.operator cimport dereference as deref
from libcpp cimport bool
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.memory cimport make_shared
cimport parpydtk2 as dtk
from mpi4py cimport libmpi as cMPI
from mpi4py cimport MPI as pyMPI


ctypedef pyMPI.Comm comm_t
ctypedef cMPI.MPI_Comm c_comm_t


cdef extern from 'mpi-compat.hpp':
    pass


import numpy as np
from ._version import __version__


__version__ = __version__
__author__ = 'Qiao Chen'
__copyright__ = 'Copyright 2018, Qiao Chen'


cdef class IMeshDB(object):
    """Interface mesh database

    ParPyDTK2 utilizes MOAB as the underlying mesh database. MOAB is an array
    based mesh library that is adapted by DTK2. With array based mesh library,
    the memory usage and computational cost are lower than typical pointer
    based data structure. The mesh concept in this work is simple since only
    meshless methods are ultilized, the only additional attribute one needs
    is the `global IDs/handles`, which are used by both MOAB and DTK2. For most
    applications, the global IDs can be computed offline.

    One thing is not directly supported by IMeshDB is I/O. However, since this
    is a Python module and only points clouds are needed, one can easily uses
    a tool (e.g. `meshio <https://github.com/nschloe/meshio>`_) to load the
    mesh.

    Attributes
    ----------
    comm : MPI.Comm
        MPI communicator
    ranks : int
        size of comm
    rank : int
        rank of comm
    size : int
        point cloud size, i.e. number of vertices
    bbox : np.ndarray
        local bounding box array of shape (2,3)
    gbbox : np.ndarray
        global bounding box array of shape (2,3)
    """

    def __init__(self, comm=None):
        """Constructor
        
        Parameters
        ----------
        comm : MPI.Comm (optional)
            if no communicator or ``None`` is passed in, then ``MPI_COMM_WORLD``
            will be used

        Examples
        --------
        >>> # implicit communciator
        >>> import parpydtk2 as dtk
        >>> mdb = dtk.IMeshDB()

        >>> # explicit communicator
        >>> from mpi4py import MPI
        >>> import parpydtk2 as dtk
        >>> mdb = dtk.IMeshDB(MPI.COMM_WOLRD)
        """
        pass

    def __cinit__(self, comm_t comm=None):
        cdef c_comm_t comm_
        if comm is None:
            comm_ = cMPI.MPI_COMM_WORLD
        else:
            comm_ = comm.ob_mpi
        self.mdb = make_shared[dtk.IMeshDB](comm_)

    def __dealloc__(self):
        self.mdb.reset()

    @property
    def comm(self):
        """MPI.Comm: communicator"""
        cdef comm_t comm = pyMPI.Comm()
        comm.ob_mpi = self.mdb.get().comm()
        return comm

    @property
    def ranks(self):
        """int: Get the communicator size"""
        return self.mdb.get().ranks()

    @property
    def rank(self):
        """int: get the rank"""
        return self.mdb.get().rank()

    def begin_create(self):
        """Begin to create/manupilate the mesh

        This function must be called in order to let the mesh databse be
        aware that you will create meshes.

        See Also
        --------
        :func:`finish_create` : finish creating mesh
        """
        self.mdb.get().begin_create()

    def create_vertices(self, double[:, ::1] coords not None):
        """Create a set of coordinates

        .. note:: The ``coords`` must be C-ordering with ndim=2!

        Parameters
        ----------
        coords : np.ndarray
            nx3 coordinates in double precision

        See Also
        --------
        :func:`assign_gids` : assign global IDs
        :func:`extract_vertices`: extract vertex coordinates

        Examples
        --------
        >>> from parpydtk2 import *
        >>> import numpy as np
        >>> mdb1 = IMeshDB()
        >>> mdb1.begin_create()
        >>> verts = np.zeros((2,3))  # two nodes
        >>> verts[1][0] = 1.0
        >>> mdb1.create_vertices(verts)
        """
        assert coords.shape[1] == 3
        self.mdb.get().create_vertices(
            <int> coords.shape[0],
            <const double *> &coords[0, 0]
        )

    def extract_vertices(self):
        """Extract coordinate
        
        .. warning::
        
            This function should be called once you have finished
            :func:`create_vertices`.

        Returns
        -------
        np.ndarray:
            (nx3) array that stores the coordinate values

        See Also
        --------
        :attr:`size` : get the mesh size
        """
        cdef cnp.ndarray[double, ndim=2] co = np.empty((self.mdb.get().size(), 3))
        self.mdb.get().extract_vertices(<double *> co.data)
        return co

    def assign_gids(self, int[::1] gids):
        """Assign global IDs

        Internally, both DTK and MOAB use so-called global IDs/handles
        communications. Each node has its own local IDs/handles and a unique
        global ID.

        Parameters
        ----------
        gids : np.ndarray
            global IDs

        See Also
        --------
        :func:`create_vertices` : create vertices
        :func:`extract_gids` : extract global IDs
        """
        self.mdb.get().assign_gids(<int> gids.size, <const int *> &gids[0])

    def extract_gids(self):
        """Extract global IDs/handles
        
        .. warning::
        
            This function should be called once you have finished
            :func:`assign_gids`.

        Returns
        -------
        np.ndarray:
            array of :attr:`size` that stores the integer IDs
        """
        cdef cnp.ndarray[int, ndim=1] gids = \
            np.empty(self.mdb.get().size(), dtype='intc')
        self.mdb.get().extract_gids(<int *> gids.data)
        return gids

    def finish_create(self, trivial_gid=True):
        """finish mesh creation

        This method finalizes the interface mesh database by communicating
        the bounding boxes and empty partitions. Also, setting up the DTK
        managers happens here.

        .. warning::
        
            You must call this function once you have done with manupilating
            the mesh, i.e. vertices and global IDs.

        Parameters
        ----------
        trivial_gid : bool
            `True` if we use MOAB trivial global ID computation

        Notes
        -----
        By ``trivial_gid``, it means simply assigning the global IDs based on
        the size of the mesh. This is useful in serial settings or transferring
        solutions from a serial solver to a partitioned one.
        """
        cdef bool tg = <bool> 1 if trivial_gid else <bool> 0
        self.mdb.get().finish_create(tg)

    @property
    def size(self):
        """int: Get the size of a set"""
        return self.mdb.get().size()

    def empty(self):
        """Check if this is an empty partition"""
        return self.mdb.get().empty()

    def has_empty(self):
        """Check if an empty partition exists"""
        return self.mdb.get().has_empty()

    @property
    def bbox(self):
        """np.ndarray: local bounding box

        The bounding box is stored simply in a 2x3 array, where the first row
        stores the maximum bounds while minimum bounds for the second row.

        .. warning::

            Bounding box is valid only after :func:`finish_create`.
        
        See Also
        --------
        :attr:`gbbox`: global bounding box
        """
        cdef cnp.ndarray[double, ndim=2, mode='c'] box = \
            np.empty(shape=(2, 3), dtype='double')
        self.mdb.get().get_bbox(<double *> box.data)
        return box

    @property
    def gbbox(self):
        """np.ndarray: global bounding box
        
        The bounding box is stored simply in a 2x3 array, where the first row
        stores the maximum bounds while minimum bounds for the second row.

        .. warning::

            Bounding box is valid only after :func:`finish_create`.

        See Also
        --------
        :attr:`bbox`: local bounding box
        """
        cdef cnp.ndarray[double, ndim=2, mode='c'] box = \
            np.empty(shape=(2, 3), dtype='double')
        self.mdb.get().get_gbbox(<double *> box.data)
        return box

    def create_field(self, str field_name, int dim=1):
        """Create a data field for solution transfer

        This is the core function to register a field so that you can then
        transfer its values to other domains. The ``dim`` parameter determines
        the data type of the field. By default, it's 1, i.e. scalar fields.
        For each node, a tensor of (1x``dim``) can be registered. For instance,
        to transfer forces and displacements in FSI applications, ``dim`` is
        3 (for 3D problems).

        Parameters
        ----------
        field_name : str
            name of the field
        dim : int
            dimension of the field, i.e. scalar, vector, tensor

        Examples
        --------
        >>> from parpydtk2 import *
        >>> mdb1 = IMeshDB()
        >>> mdb1.begin_create()
        >>> mdb1.create_field('heat flux')
        """
        cdef std_string fn = <std_string> field_name.encode('UTF-8')
        self.mdb.get().create_field(fn, dim)

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
        return self.mdb.get().has_field(fn)

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
        return self.mdb.get().field_dim(fn)

    def assign_field(self, str field_name, cnp.ndarray values not None):
        """Assign values to a field

        .. note:: `values` size must be at least size*dim

        Parameters
        ----------
        field_name : str
            name of the field
        values : np.ndarray
            input source values

        See Also
        --------
        :func:`extract_field` : extract value from a field
        :attr:`size` : check the size of a mesh set
        """
        cdef:
            std_string fn = <std_string> field_name.encode('UTF-8')
            int dim = self.mdb.get().field_dim(fn)
            int n = self.mdb.get().size() * dim
        assert values.size >= n
        self.mdb.get().assign_field(fn, <const double *> values.data)

    def resolve_empty_partitions(self, str field_name):
        """Resolve asynchronous values on empty partitions

        ParPyDTK2 doesn't expect :func:`assign_field` should be called
        collectively. Therefore, a collective call must be made for resolving
        empty partitions.

        .. warning:: This must be called collectively even on empty partitions

        .. note:: You should call this function following assignment

        Parameters
        ----------
        field_name : str
            name of the field

        Examples
        --------
        >>> if rank == 0:
        ...     mdb.assign_field('flux', values)
        >>> mdb.resolve_empty_partitions('flux')

        Notes
        -----
        This API is not available in C++ level, therefore, one needs to
        implement this if he/she wants to use the C++ API.
        """
        cdef:
            std_string fn = <std_string> field_name.encode('UTF-8')
            int dim = self.mdb.get().field_dim(fn)
            const std_vector[int] *m2s
            int i
            cnp.ndarray[double, ndim=1] buf
        if self.mdb.get().has_empty():
            buf = np.empty(dim)
            if self.mdb.get().rank() == 0:
                m2s = &self.mdb.get()._m2s()
                self.mdb.get()._extract_1st(fn, <double *> buf.data)
                for i in range(deref(m2s).size()):
                    self.comm.Isend(buf, dest=deref(m2s)[i])
            elif self.mdb.get().empty():
                req = self.comm.Irecv(buf, source=0)
                req.Wait()
                self.mdb.get()._assign_1st(fn, <const double *> buf.data)


    def extract_field(self, str field_name,
        double[::1] buffer=None, reshape=False):
        """Extact the values from a field

        .. warning:: if `buffer` is passed in, it must be 1D

        Parameters
        ----------
        field_name : str
            name of the field
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
            int n = self.mdb.get().size()
            int dim = self.mdb.get().field_dim(fn)
            cnp.ndarray[double, ndim=1] values
        if buffer is not None:
            assert len(buffer) >= n * dim
            values = np.asarray(buffer)
        else:
            values = np.empty(n * dim, dtype='double')
        self.mdb.get().extract_field(fn, <double *> values.data)
        if reshape and dim > 1:
            return values.reshape((n, dim))
        return values
