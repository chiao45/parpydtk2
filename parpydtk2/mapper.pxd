"""Cython interface of Python modules"""

cimport parpydtk2 as dtk
from .imeshdb cimport IMeshDB


cdef class Mapper(object):
    cdef dtk.Mapper *mp
    cdef IMeshDB blue
    cdef IMeshDB green
