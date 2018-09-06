"""Cython interface of IMeshDB"""

cimport parpydtk2 as dtk
from libcpp.memory cimport shared_ptr


cdef class IMeshDB(object):
    cdef shared_ptr[dtk.IMeshDB] mdb