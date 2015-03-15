# distutils: language = c++
from cpython.array cimport array as pyarray
from pytraj.cpptraj_dict import MatrixDict, MatrixKindDict, get_key


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

    @property
    def mkind(self):
        return get_key(self.baseptr_1.Kind(), MatrixKindDict)
