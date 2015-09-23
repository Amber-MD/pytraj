# distutils: language = c++

from libcpp.string cimport string
from ..core.cpptraj_core cimport CpptrajState, _CpptrajState

cdef class DataFile:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFile()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def read_data(self, filenameIn, arglist, DatasetList datasetlist):
        return self.thisptr.ReadDataIn(filenameIn.encode(),
               ArgList(arglist).thisptr[0], datasetlist.thisptr[0])


cdef class DataFileList:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFileList()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def write_all_datafiles(self):
        self.thisptr.WriteAllDF()
