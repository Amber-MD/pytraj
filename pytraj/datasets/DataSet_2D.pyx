# distutils: language = c++
from cpython.array cimport array as pyarray
from ..cpptraj_dict import MatrixDict, MatrixKindDict, get_key


cdef class DataSet_2D (DataSet):
    def __cinit__(self):
        # since DataSet_2D inherits from DataSet, make sure two pointers pointing 
        # to the same address
        self.baseptr_1 = <_DataSet_2D*> self.baseptr0

    def __dealloc__(self):
        pass

    @property
    def n_rows(self):
        return self.baseptr_1.Nrows()

    @property
    def n_cols(self):
        return self.baseptr_1.Ncols()

    def get_element(self, int x, int y):
        return self.baseptr_1.GetElement(x, y)

    @property
    def mkind(self):
        """return matrix kind: full, half or triangle"""
        return get_key(self.baseptr_1.Kind(), MatrixKindDict)

    def allocate_2D(self, size_t x, size_t y):
        self.baseptr_1.Allocate2D(x, y)

    def allocate_half(self, size_t x):
        self.baseptr_1.AllocateHalf(x)

    def allocate_triangle(self, size_t x):
        self.baseptr_1.AllocateTriangle(x)

    def get_full_matrix(self):
        raise NotImplementedError("must over-write in subclass")
