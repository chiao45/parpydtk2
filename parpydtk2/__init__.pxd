"""ParPyDTK2 Cython interface"""

from libcpp cimport bool
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.memory cimport shared_ptr
from mpi4py cimport libmpi as cMPI


ctypedef cMPI.MPI_Comm c_comm_t


cdef extern from 'src/dtk2.hpp' namespace 'parpydtk2' nogil:

    cdef cppclass IMeshDB:
        IMeshDB(c_comm_t comm) except +
        int ranks()
        int rank()
        c_comm_t comm()
        void begin_create() except +
        void create_vset() except +
        void create_vertices(int nv, const double *coords) except +
        void extract_vertices(double *coords) except +
        void assign_gids(int nv, const int *gids) except +
        void extract_gids(int *gids) except +
        void finish_create(bool trivial_gid) except +
        int size()
        bool empty()
        bool has_empty()
        const std_vector[int] &_m2s()
        void get_bbox(double *v) except +
        void get_gbbox(double *v) except +
        void create_field(const std_string &field_name, int dim) except +
        bool has_field(const std_string &field_name)
        int field_dim(const std_string &field_name) except +
        void assign_field(
            const std_string &field_name,
            const double *values
        ) except +
        void _assign_1st(
            const std_string &field_name,
            const double *values
        ) except +
        void extract_field(
            const std_string &field_name,
            double *values
        ) except +
        void _extract_1st(
            const std_string &field_name,
            double *values
        ) except +


    cdef cppclass Mapper:
        Mapper(
            shared_ptr[IMeshDB] B,
            shared_ptr[IMeshDB] G,
            const std_string &version,
            const std_string &date,
            bool profiling
        ) except +
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
        shared_ptr[IMeshDB] blue_mesh()
        shared_ptr[IMeshDB] green_mesh()
        void begin_initialization() except +
        void register_coupling_fields(
            const std_string &bf,
            const std_string &gf,
            bool direct
        ) except +
        bool has_coupling_fields(
            const std_string &bf,
            const std_string &gf,
            bool direct
        )
        void end_initialization() except +
        void begin_transfer()
        void transfer_data(
            const std_string &bf,
            const std_string &gf,
            bool direct
        ) except +
        void end_transfer()
