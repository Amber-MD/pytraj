# distutils: language = c++
from cpython.array cimport array as pyarray
from ..cpptraj_dict import MatrixDict, MatrixKindDict, get_key

cdef class DataSet_MatrixFlt (DataSet_2D):
    def __cinit__(self):
        self.thisptr = new _DataSet_MatrixFlt()
        self.baseptr_1 = <_DataSet_2D*> self.thisptr
        self.baseptr0 = <_DataSet*> self.thisptr

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def __getitem__(self, idx):
        return self.data[idx]

    def alloc(self):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = _DataSet_MatrixFlt.Alloc()
        return dset

    #@property
    #def mass(self):
    #    return self.thisptr.Mass()

    def get_full_matrix(self):
        """return python array with length = n_rows*n_cols"""
        cdef int nr = self.n_rows
        cdef int nc = self.n_cols 
        cdef int i, j
        cdef pyarray arr0 = pyarray('f', [])

        for i in range(nr):
            for j in range(nc):
                arr0.append(self.baseptr_1.GetElement(i, j))
        return arr0

    @property
    def data(self):
        """return 1D python array of matrix' data"""
        return self.get_full_matrix()
