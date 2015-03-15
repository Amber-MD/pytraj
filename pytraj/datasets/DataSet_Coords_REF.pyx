# distutils: language = c++


cdef class DataSet_Coords_REF (DataSet_Coords):
    def __cinit__(self):
        self.thisptr = new _DataSet_Coords_REF()
        self.baseptr0 = <_DataSet*> self.thisptr
        self.baseptr_2 = <_DataSet_Coords*> self.thisptr
        self.baseptr_1 = <_DataSet_1D*> self.thisptr

        # let python frees memory
        self.py_free_mem = True

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    @classmethod
    def alloc(self):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = _DataSet_Coords_REF.Alloc()
        return dset

    @property
    def size(self):
        return self.thisptr.Size()

    def get_frame(self):
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.thisptr.RefFrame()
        return frame
