#!python
#cython: language_level=3
#cython: boundscheck=False
#cython: embedsignature=True
#cython: wraparound=False

"""Interface solution transfer/mapper

This module wraps DTK2 meshless solution transfer operators, which can provide
convenient data remapping services.

.. moduleauthor:: Qiao Chen <benechiao@gmail.com>
"""

cimport numpy as cnp
from libcpp cimport bool
from libcpp.string cimport string as std_string
cimport parpydtk2 as dtk
from mpi4py cimport MPI as pyMPI
from .imeshdb cimport IMeshDB


ctypedef pyMPI.Comm comm_t


cdef extern from 'mpi-compat.hpp':
    pass


import numpy as np
import datetime
from ._version import __version__


__version__ = __version__
__author__ = 'Qiao Chen'
__copyright__ = 'Copyright 2018, Qiao Chen'


cdef class Mapper(object):
    """DTK2 wrapper

    The meshless methods in DTK2, including `modified moving least square`,
    `spline interpolation` and `nearest node projection` methods are wrapped
    within this class. The most advanced method is the MMLS fitting, which
    is my personal recommendataion.

    Attributes
    ----------
    blue_mesh : :py:class:`~parpydtk2.IMeshDB`
        blue mesh participant
    green_mesh : :py:class:`~parpydtk2.IMeshDB`
        green mesh participant
    comm : MPI.Comm
        MPI communicator
    ranks : int
        size of comm
    rank : int
        rank of comm
    dimension : int
        spatial dimension
    method : int
        method flag, either MMLS, SPLINE, or N2N
    basis : int
        flag of basis function for weighting schemes used by MMLS and SPLINE
    knn_b : int
        k-nearest neighborhood for searching on blue_mesh
    knn_g : int
        k-nearest neighborhood for searching on green_mesh
    radius_b : float
        radius used for searching on blue_mesh
    radius_g : float
        radius used for searching on green_mesh
    """

    def __init__(self, blue, green, profiling=True):
        """Constructor

        Parameters
        ----------
        blue_mesh : :py:class:`~parpydtk2.IMeshDB`
            blue mesh participant
        green_mesh : :py:class:`~parpydtk2.IMeshDB`
            green mesh participant
        profiling : bool (optional)
            whether or not do timing report, default is ``True``.

        Examples
        --------
        >>> from parpydtk2 import *
        >>> blue, green = IMeshDB(), IMeshDB()
        >>> # initialize blue and green
        >>> mapper = Mapper(blue=blue,green=green)
        >>> # do work with mapper
        """
        pass
    
    def __cinit__(self, IMeshDB blue, IMeshDB green, profiling=True):
        cdef:
            std_string version = __version__.encode('UTF-8')
            std_string date = \
                datetime.datetime.now().strftime('%b %d %Y %H:%M:%S').encode('UTF-8')
            bool prof = <bool> 1 if profiling else <bool> 0
        self.mp = new dtk.Mapper(blue.mdb, green.mdb, version, date, prof)
        self.blue = blue
        self.green = green
    
    def __dealloc__(self):
        del self.mp

    @property
    def blue_mesh(self):
        """:py:class:`~parpydtk2.IMeshDB`: blue mesh"""
        return self.blue

    @property
    def green_mesh(self):
        """:py:class:`~parpydtk2.IMeshDB`: green mesh"""
        return self.green

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
        cdef comm_t comm = pyMPI.Comm()
        comm.ob_mpi = self.mp.comm()
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
            ``True`` if (bf,gf) exists in the ``direct`` direction
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
        """Transfer (bf, gf) in the ``direct`` direction

        Parameters
        ----------
        bf : str
            blue mesh field name
        gf : str
            green mesh field name
        direct : bool
            ``True`` for blue->green, ``False`` for the opposite

        See Also
        --------
        :func:`register_coupling_fields` : register coupled fields
        """
        cdef:
            std_string bf_ = <std_string> bf.encode('UTF-8')
            std_string gf_ = <std_string> gf.encode('UTF-8')
            bool dr = <bool> 1 if direct else <bool> 0
        # resolve empty partitions
        if dr:
            src = self.blue
            tag_src = bf
        else:
            src = self.green
            tag_src = gf
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