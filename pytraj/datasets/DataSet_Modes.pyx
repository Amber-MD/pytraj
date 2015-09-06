# distutils: language = c++


cdef class DataSet_Modes (DataSet):
    def __cinit__(self):
        self.thisptr = new _DataSet_Modes()
        self.baseptr0 = <_DataSet*> self.thisptr

    def __dealloc__(self):
        del self.thisptr

    def alloc(self):
        '''return a memoryview as DataSet instane'''
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.Alloc()
        return dset


    def nmodes(self):
        return self.thisptr.Nmodes()

    def vector_size(self):
        return self.thisptr.VectorSize()

    def is_reduced(self):
        return self.thisptr.IsReduced()

