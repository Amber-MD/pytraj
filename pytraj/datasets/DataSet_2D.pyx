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
    def data(self):
        cdef int i
        cdef pyarray arr = pyarray('d', [])

        if self.baseptr_1.MatrixArray():
            for i in range(self.n_cols * self.n_rows):
                arr.append(self.baseptr_1.MatrixArray()[i])
            return arr
        else:
            raise ValueError("Can not get MatrixArray")

    @property
    def type(self):
        return get_key(self.baseptr_1.m2dType(), MatrixDict)
