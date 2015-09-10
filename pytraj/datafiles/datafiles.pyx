# distutils: language = c++

from pytraj.cpptraj_dict import DataFormatDict, get_key
from pytraj.externals.six import string_types
from cython.operator cimport dereference as deref


cdef class DataFile:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFile()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

cdef class DataFileList:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFileList()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr
