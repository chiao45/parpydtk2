"""DTK2Mapper cython file"""

from libcpp.string cimport string as std_string
from libcpp cimport bool
from pymoab cimport moab
from pymoab.core cimport Core


cdef extern from 'DtkWrap.hpp' namespace 'pydtk2':

    cdef cppclass Dtk2MoabManager:

        Dtk2MoabManager(moab.Interface *meshdb, unsigned long mesh_set) except +

        void register_tag(const std_string &tag) except +
        bool has_tag(const std_string &tag)


    cdef cppclass Dtk2Mapper:

        Dtk2Mapper(
            Dtk2MoabManager *source,
            Dtk2MoabManager *target,
            double r,
            int basis,
            int order,
            int dim,
            bool track_missing) except +

        void set_basis(int basis) except +
        void set_spacial_dimension(int dim) except +
        void set_order(int order) except +
        void set_rbf_radius(double r) except +
        void set_track_missing_flag(bool flag) except +

        void register_coupled_tags(
            const std_string &stag, const std_string &ttag) except +
        void apply(
            const std_string &stag,
            const std_string &ttag) except +
        bool has_coupled_tags(const std_string &stag, const std_string &ttag)


cdef class DTK2MoabManager:

    cdef Dtk2MoabManager *_inst

    def __cinit__(self, Core meshdb, mesh_set=0):
        """Constructor"""
        cdef unsigned int mset = <unsigned int> mesh_set
        self._inst = new Dtk2MoabManager(<moab.Interface *> meshdb.inst, mset)

    def __del__(self):
        """Destructor"""
        del self._inst

    def register_tag(self, str tag):
        """Register a tag to mapper,

        The tag string must be already created in MOAB core
        """
        cdef std_string t = <std_string> tag.encode('UTF-8')
        self._inst.register_tag(t)

    def has_tag(self, str tag):
        """Check if a tag is already registered"""
        cdef std_string t = <std_string> tag.encode('UTF-8')
        return self._inst.has_tag(t)

    # pure python
    def register_tags(self, tags):
        """Register a list of tags"""
        for tag in tags:
            self.register_tag(tag)


cdef class DTK2Mapper:

    cdef Dtk2Mapper *_inst

    def __cinit__(self,
                  DTK2MoabManager source,
                  DTK2MoabManager target,
                  double r,
                  basis_order=None,
                  dim=3,
                  track_missing=True):
        """Constructor"""
        cdef double radius = <double> r
        cdef bool tm = <bool> track_missing
        if basis_order is None:
            basis = 0
            order = 4
        else:
            basis = basis_order[0]
            order = basis_order[1]
        if [basis, order] not in [[0, 0], [0, 2], [0, 4], [0, 6],
                                  [1, 0], [1, 2], [1, 4], [2, 3]]:
            raise ValueError('Unsupported BASIS/ORDER.')
        self._inst = new Dtk2Mapper(
            <Dtk2MoabManager *> source._inst,
            <Dtk2MoabManager *> target._inst,
            radius,
            <int> basis,
            <int> order,
            <int> dim,
            tm)

    def __del__(self):
        """Destructor"""
        del self._inst

    def set_spacial_dimension(self, int dim):
        """Set spatial dimension"""
        assert dim > 0 and dim < 4, 'Invalid dimension.'
        self._inst.set_spacial_dimension(<int> dim)

    def set_basis_order(self, basis_order):
        """Set basis and order"""
        assert basis_order in [
            [0, 0], [0, 2], [0, 4], [0, 6],
            [1, 0], [1, 2], [1, 4], [2, 3]], 'Unknown BASIS/ORDER.'
        self._inst.set_basis(<int> basis_order[0])
        self._inst.set_order(<int> basis_order[1])

    def set_rbf_radius(self, double r):
        """Set the radias basis functions radius"""
        self._inst.set_rbf_radius(<double> r)

    def set_track_missing_flag(self, bool flag):
        """Set the track missing elements flag"""
        self._inst.set_track_missing_flag(<bool> flag)

    def register_coupled_tags(self, str stag, str ttag):
        """Register a coupled tags into this mapper

        The tags must be already registered in managers
        for source and target respectively.
        """
        cdef std_string st = <std_string> stag.encode('UTF-8')
        cdef std_string tt = <std_string> ttag.encode('UTF-8')
        self._inst.register_coupled_tags(st, tt)

    def apply(self, str stag, str ttag):
        """Apply the solution transfer"""
        cdef std_string st = <std_string> stag.encode('UTF-8')
        cdef std_string tt = <std_string> ttag.encode('UTF-8')
        self._inst.apply(st, tt)

    def has_coupled_tags(self, str stag, str ttag):
        """Check if a coupled tags of <stag, ttag> is already regiestered"""
        cdef std_string st = <std_string> stag.encode('UTF-8')
        cdef std_string tt = <std_string> ttag.encode('UTF-8')
        return self._inst.has_coupled_tags(st, tt)
