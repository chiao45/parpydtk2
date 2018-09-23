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


_ATTR_PARS = {
    'method',
    'basis',
    'dimension',
    'radius_b',
    'radius_g',
    'leaf_b',
    'leaf_g',
    'rho'
}


cdef class Mapper(object):
    """DTK2 wrapper

    The meshless methods in DTK2, including `modified moving least square`,
    `spline interpolation` and `nearest node projection` methods are wrapped
    within this class.

    In addition, if you build DTK2 from UNIFEM/CHIAO45 forked verions, then
    you can use the `adaptive weighted least square` fitting method.

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
    radius_b : float
        radius used for searching on blue_mesh
    radius_g : float
        radius used for searching on green_mesh
    leaf_b : int
        kd-tree leaf size of blue_mesh
    leaf_g : int
        kd-tree leaf size of green_mesh
    rho : double
        local Vandermonde system row scaling factor
    """

    @staticmethod
    def parameter_keys():
        from copy import deepcopy
        keys = deepcopy(_ATTR_PARS)
        keys.add('resolve_disc')
        keys.add('_ind_file')
        return keys

    @staticmethod
    def is_unifem_backend():
        """Check if the underlying DTK2 library is unifem forked version

        Returns
        -------
        bool
            ``True`` if the user's DTK2 is from unifem forked repo
        """
        return dtk.Mapper.is_unifem_backend()

    def __init__(self, blue, green, profiling=True, verbose=True, **kwargs):
        """Constructor

        Parameters
        ----------
        blue_mesh : :py:class:`~parpydtk2.IMeshDB`
            blue mesh participant
        green_mesh : :py:class:`~parpydtk2.IMeshDB`
            green mesh participant
        profiling : bool (optional)
            whether or not do timing report, default is ``True``.
        verbose : bool (optional)
            whether or not verbose printing, default is ``True``.
        stat_file : str (optional)
            file for storing profiling information

        Examples
        --------
        >>> from parpydtk2 import *
        >>> blue, green = IMeshDB(), IMeshDB()
        >>> # initialize blue and green
        >>> mapper = Mapper(blue=blue,green=green)
        >>> # do work with mapper

        ParPyDTK2 Profiling File
        ------------------------

        By default, if ``profiling`` is required, a mapper will create an
        ASCII file named ``mapper_[version].stat``, where ``mapper_[version]``
        can be specified by the user with parameter ``stat_file``. If the
        program runs in parallel, then a suffix of ``_[rank]`` will be added
        to the filename, i.e. ``mapper_[version]_[rank].stat`` or
        ``[stat_file]_[rank].stat``.
        """
        pass

    def __cinit__(self, IMeshDB blue, IMeshDB green, *,
        profiling=True, str stat_file='', verbose=True, **kwargs):
        cdef:
            std_string version = __version__.encode('UTF-8')
            std_string date = \
                datetime.datetime.now().strftime('%b %d %Y %H:%M:%S').encode('UTF-8')
            bool prof = <bool> 1 if profiling else <bool> 0
            bool verb = <bool> 1 if verbose else <bool> 0
            std_string sf = stat_file.encode('UTF-8')
        self.mp = new dtk.Mapper(blue.mdb, green.mdb, version, date, prof, sf, verb)
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
        elif m == 3:
            if not Mapper.is_unifem_backend():
                import warnings
                warnings.warn(
                    'The underlying DTK2 is not from UNIFEM or CHIAO45, AWLS will reduce to MMLS',
                    RuntimeWarning
                )
            self.mp.use_awls()
        elif m == 4:
            self.mp.use_n2n(<bool> 1)
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
        import warnings
        warnings.warn(
            'This method will be deprecated, direct assign method with N2N_MATCH',
            DeprecationWarning
        )
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

    @property
    def leaf_b(self):
        """int: get the leaf size of blue mesh for kd-tree
        
        .. warning::

            This attribute only works when the underlying DTK2 is installed
            from UNIFEM or CHIAO45 forked versions.
        
        See Also
        --------
        :attr:`leaf_g` : green leaf size of kd-tree
        """
        return self.mp.leaf_b()
    
    @leaf_b.setter
    def leaf_b(self, int size):
        if not Mapper.is_unifem_backend():
            import warnings
            warnings.warn(
                'Leaf size is only tunable while the underlying DTK2 is from UNIFEM or CHIAO45 forked versions',
                RuntimeWarning
            )
        self.mp.set_leaf_b(size)
        if size < 5 or size > 50:
            import warnings
            warnings.warn(
                'Too large/small leaf size:{}'.format(size),
                RuntimeWarning
            )
    
    @property
    def leaf_g(self):
        """int: get the leaf size of the green mesh for kd-tree
        
        .. warning::

            This attribute only works when the underlying DTK2 is installed
            from UNIFEM or CHIAO45 forked versions.

        See Also
        --------
        :attr:`leaf_b` : blue leaf size of kd-tree
        """
        return self.mp.leaf_g()
    
    @leaf_g.setter
    def leaf_g(self, int size):
        if not Mapper.is_unifem_backend():
            import warnings
            warnings.warn(
                'Leaf size is only tunable while the underlying DTK2 is from UNIFEM or CHIAO45 forked versions',
                RuntimeWarning
            )
        self.mp.set_leaf_g(size)
        if size < 5 or size > 50:
            import warnings
            warnings.warn(
                'Too large/small leaf size:{}'.format(size),
                RuntimeWarning
            )
    
    def _set_ind_file(self, str filename):
        # WARNING! should be used in serial for tuning parameter!
        cdef std_string fn = filename.encode('UTF-8')
        import warnings
        if not Mapper.is_unifem_backend():
            warnings.warn(
                'The underlying DTK2 installation is not from UNIFEM!',
                RuntimeWarning
            )
        self.mp._set_ind_file(fn)
    
    def _wipe_ind_file(self):
        self.mp._wipe_ind_file()
    
    @property
    def rho(self):
        "float: local scaling factor"
        return self.mp.rho()

    @rho.setter
    def rho(self, double v):
        self.mp.set_rho(v)

    def enable_unifem_mmls_auto_conf(self, *, ref_r_b=None, ref_r_g=None,
        dim=None, verbose=False, **kwargs):
        r"""Automatically set up radius parameter for MMLS

        .. warning::

            This method should be used only when the underlying DTK2
            installation is from UNIFEM or CHIAO45 forked versions. Otherwise,
            you should **always** manully configure the radius parameters.

        The following strategy is performed:

        .. math::

            r = \max(\alpha h_b, \frac{\beta h_b}{N^{1/d}}, r_u)

        where :math:`h_b` is the the the maximum edge length of the global
        bounding box (:py:attr:`~parpydtk2.IMeshDB.gbbox`); :math:`\alpha` is
        some ratio, say 0.1 (10%); :math:`\beta` is the scaling factor the
        the estimated mesh size, which is given by
        :math:`\frac{h_b}{N^{1/(d-1)}}` and :math:`N` is the global mesh size
        (:py:attr:`~parpydtk2.IMeshDB.gsize`); :math:`d` is the spatial
        dimension; the last parameter :math:`r_u` is provided by the user thus
        optional.

        For UNIFEM/CHIAO45 DTK2, this parameter should be relatively large,
        because the final points in the Vandermonde system is determined by
        the column size, so that larger radius means that the system has
        larger candidate pool.

        .. warning::

            Regarding the spatial dimension, since this packages is mainly for
            interface/surface coupling thus the actual dimension is assumed to
            be one less than the spatial dimension, i.e. surface topological
            dimension. If this is not the case, the user needs to explicit pass
            in the dimension to override this default behavior.

        Parameters
        ----------
        ref_r_b : float (optional)
            reference user-specified blue radius, i.e. :math:`r_u` for blue
        ref_r_g : float (optional)
            reference user-specified green radius, i.e. :math:`r_u` for green
        dim : int (optional)
            topological dimension, default is the surface dimension
        verbose : bool (optional)
            print verbose information/warning messages, default is ``False``
        alpha : float
            the :math:`\alpha` parameter
        beta : float
            the :math:`\beta` parameter

        See Also
        --------
        :func:`is_unifem_backend` : check backend installation of DTK2
        """
        import warnings
        import sys
        if not Mapper.is_unifem_backend():
            warnings.warn(
                'The underlying DTK2 installation is not from UNIFEM!',
                RuntimeWarning
            )
        alpha = kwargs.pop('alpha', 0.1)
        if alpha <= 0.0:
            if verbose:
                warnings.warn('Invalid alpha:{}'.format(alpha), RuntimeWarning)
            alpha = 0.1
        beta = kwargs.pop('beta', 5.0)
        if beta <= 0.0:
            if verbose:
                warnings.warn('Invalid beta:{}'.format(beta), RuntimeWarning)
            beta = 5.0
        # first check the method
        self.method = 3  # set to awls
        if dim is None:
            dim = self.dimension - 1  # assume surface
        elif int(dim) < 0 or int(dim) > 3:
            if verbose:
                warnings.warn('Invalid dimension %r' % dim, RuntimeWarning)
            dim = self.dimension - 1  # assume surface
        else:
            dim = int(dim)
        if verbose:
            print('Mapper dimension is: {}.'.format(dim))

        # blue radius
        b_gsize = self.blue_mesh.gsize
        b_gbox = self.blue_mesh.gbbox
        b_h = -1.0
        for i in range(self.dimension):
            b_h = max(b_h, abs(b_gbox[0][i]-b_gbox[1][i]))
        b_r1 = b_h * alpha
        b_r2 = beta * b_h / np.power(b_gsize, 1./dim)
        if ref_r_b is None:
            ref_r_b = 0.0
        elif float(ref_r_b) <= 0.0:
            if verbose:
                warnings.warn(
                    'Invalid blue user-provided radius %r' % ref_r_b,
                    RuntimeWarning
                )
            ref_r_b = 0.0
        else:
            ref_r_b = float(ref_r_b)
        r_b = max(max(b_r1, b_r2), ref_r_b)
        if verbose:
            print('Blue radius-conf: r_box={}, r_h={}, r_user={}'.format(b_r1, b_r2, ref_r_b))
            print('Blue final radius={}'.format(r_b))
        self.radius_b = r_b

        # green radius
        g_gsize = self.green_mesh.gsize
        g_gbox = self.green_mesh.gbbox
        g_h = -1.0
        for i in range(self.dimension):
            g_h = max(g_h, abs(g_gbox[0][i]-g_gbox[1][i]))
        g_r1 = g_h * alpha
        g_r2 = beta * g_h / np.power(g_gsize, 1./dim)
        if ref_r_g is None:
            ref_r_g = 0.0
        elif float(ref_r_g) <= 0.0:
            if verbose:
                warnings.warn(
                    'Invalid green user-provided radius %r' % ref_r_g,
                    RuntimeWarning
                )
            ref_r_g = 0.0
        else:
            ref_r_g = float(ref_r_g)
        r_g = max(max(g_r1, g_r2), ref_r_g)
        if verbose:
            print('Green radius-conf: r_box={}, r_h={}, r_user={}'.format(g_r1, g_r2, ref_r_g))
            print('Green final radius={}'.format(r_g))
        self.radius_g = r_g
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        except:
            pass

    def awls_conf(self, **kwargs):
        r"""Configuration of Adaptive Weighted Least Square Fitting

        Parameters
        ----------
        basis : int (optional)
            basis weighting scheme, default is WENDLAND21
        rho : float (optional)
            number of rows in the local Vandermonde system, i.e. rho*col
        ref_r_b : float (optional)
            reference user-specified blue radius, i.e. :math:`r_u` for blue
        ref_r_g : float (optional)
            reference user-specified green radius, i.e. :math:`r_u` for green
        dim : int (optional)
            topological dimension, default is the surface dimension
        verbose : bool (optional)
            print verbose information/warning messages, default is ``False``
        alpha : float
            the :math:`\alpha` parameter
        beta : float
            the :math:`\beta` parameter
        _ind_file : str (optional)
            indicator result file, used in unifem/chiao45 dtk

        Notes
        -----

        ``_ind_file`` is for internal use to fine tune the parameter ``sigma``.
        It will not be enabled in release build.

        See Also
        --------
        :func:`enable_unifem_mmls_auto_conf` : configure radius for parallel
        mesh rendezvous
        """
        # filter the kwargs and enable mmls configuration
        self.enable_unifem_mmls_auto_conf(
            ref_r_b=kwargs.pop('ref_r_b', None),
            ref_r_g=kwargs.pop('ref_r_g', None),
            dim=kwargs.pop('dim', None),
            verbose=kwargs.pop('verbose', False),
            **kwargs
        )
        self.basis = kwargs.pop('basis', 7)
        self.method = 3
        f = kwargs.pop('_ind_file', None)
        if f:
            self._set_ind_file(f)

    def begin_initialization(self, **kwargs):
        """Initialization starter

        This is a must-call function in order to indicate mapper that you are
        about to initialize/register coupling fields

        See Also
        --------
        :func:`register_coupling_fields` : register coupled fields
        :func:`end_initialization` : finish initialization
        """
        self.mp.begin_initialization()

    def register_coupling_fields(self, str bf, str gf, bool direct, **kwargs):
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

    def end_initialization(self, **kwargs):
        """Initialization closer

        This is a must-call function in order to tell the mapper we are ready

        See Also
        --------
        :func:`begin_initialization` : initialization starter
        """
        self.mp.end_initialization()

    def begin_transfer(self, **kwargs):
        """Transfer starter

        This is a must-call function to inidate the beginning of a transferring
        block

        See Also
        --------
        :func:`end_transfer` : transfer closer
        :func:`tranfer_data` : transfer a coupled data fields
        """
        self.mp.begin_transfer()

    def transfer_data(self, str bf, str gf, bool direct, **kwargs):
        """Transfer (bf, gf) in the ``direct`` direction

        .. note::

            Parameters ``resolve_disc`` and ``sigma`` only work with
            UNIFEM/CHIAO45 DTK and AWLS method.

        Parameters
        ----------
        bf : str
            blue mesh field name
        gf : str
            green mesh field name
        direct : bool
            ``True`` for blue->green, ``False`` for the opposite
        resolve_disc : bool (optional)
            ``True`` if do post-processing for resolving non-smooth functions
        sigma : float (optional)
            threshold used in smoothness indicator

        See Also
        --------
        :func:`register_coupling_fields` : register coupled fields
        """
        cdef:
            std_string bf_ = <std_string> bf.encode('UTF-8')
            std_string gf_ = <std_string> gf.encode('UTF-8')
            bool dr = <bool> 1 if direct else <bool> 0
            bool disc = <bool> 0
            double sgm = <double> kwargs.get('sigma', -1.0)
        if kwargs.get('resolve_disc', False):
            disc = <bool> 1
        self.mp.transfer_data(bf_, gf_, dr, disc, sgm)

    def end_transfer(self, **kwargs):
        """Transfer closer

        This is a must-call function to indicate we have finished a sequence of
        transferring requests

        See Also
        --------
        :func:`begin_transfer` : transfer starter
        """
        self.mp.end_transfer()

    def __getitem__(self, str key):
        if key in _ATTR_PARS:
            return getattr(self, key)
        raise KeyError('{} is not a valid key'.format(key))

    def __setitem__(self, str key, v):
        if key in _ATTR_PARS:
            setattr(self, key, v)
        elif key == '_ind_file':
            if v:
                self._set_ind_file(v)
            else:
                self._wipe_ind_file()
        else:
            raise KeyError('{} is not a valid key'.format(key))
