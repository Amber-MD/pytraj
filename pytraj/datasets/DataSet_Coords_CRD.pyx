# distutils: language = c++


cdef class DataSet_Coords_CRD (DataSet_Coords):
    def __cinit__(self):
        self.thisptr = new _DataSet_Coords_CRD()
        self.baseptr0 = <_DataSet*> self.thisptr
        self.baseptr_2 = <_DataSet_Coords*> self.thisptr
        self.baseptr_1 = <_DataSet_1D*> self.thisptr

    def __dealloc__(self):
        del self.thisptr

    #def DataSet_Coords_CRD(self):

    @classmethod
    def alloc(self):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = _DataSet_Coords_CRD.Alloc()
        return dset

    @property
    def size(self):
        return self.thisptr.Size()
