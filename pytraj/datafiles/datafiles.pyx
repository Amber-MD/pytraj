# distutils: language = c++

from libcpp.string cimport string
from ..core.c_core cimport CpptrajState, _CpptrajState

cdef class DataFile:
    def __cinit__(self, _own_memory=True):
        self.thisptr = new _DataFile()
        self._own_memory = _own_memory

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def read_data(self, filenameIn, arglist, DatasetList datasetlist):
        return self.thisptr.ReadDataIn(
            filenameIn.encode(),
            ArgList(arglist).thisptr[0],
            datasetlist.thisptr[0])


cdef class DataFileList:
    def __cinit__(self, _own_memory=True):
        self.thisptr = new _DataFileList()
        self._own_memory = _own_memory

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def write_all_datafiles(self):
        self.thisptr.WriteAllDF()
