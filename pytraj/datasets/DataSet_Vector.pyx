# distutils: language = c++

from cython.view cimport array as cyarray


cdef class DataSet_Vector (DataSet_1D):
    def __cinit__(self):
        self.py_free_mem = True
        self.thisptr = new _DataSet_Vector()
        self.baseptr0 = <_DataSet*> self.thisptr
        self.baseptr_1= <_DataSet_1D*> self.thisptr

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def alloc(self):
        cdef DataSet d0 = DataSet()
        d0.baseptr0 = self.thisptr.Alloc()
        return d0

    def __getitem__(self, idx):
        return self.data[idx]

    def __iter__(self):
        for i in range (self.size):
            yield self[i]

    @property
    def data(self):
        """return a list of Vec3"""
        cdef Vec3 vec
        cdef cyarray cya = cyarray(shape=(self.size, 3), itemsize=sizeof(double), format="d")
        cdef int idx

        for idx in range(self.size):
            vec = Vec3()
            vec.thisptr[0] = self.thisptr.index_opr(idx)
            cya[idx, :] = vec.buffer1d[:]
        return cya
