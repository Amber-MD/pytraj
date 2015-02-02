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
        # don't __dealloc__ here since this is abstract class
        #del self.baseptr_1

    #def DataSet_2D(self,DataSet:

    #def virtual int Allocate2D(self,size_t, size_t)= 0 :

    #def virtual int AllocateHalf(self,size_t)= 0 :

    #def virtual int AllocateTriangle(self,size_t)= 0 :

    #def virtual void Write2D(self,CpptrajFile, int, int) = 0 :

    #def virtual double GetElement(self,size_t, size_t) = 0 :

    @property
    def n_rows(self):
        return self.baseptr_1.Nrows()

    @property
    def n_cols(self):
        return self.baseptr_1.Ncols()

    #def virtual double * MatrixArray(self):
    @property
    def data(self):
        cdef int i
        cdef pyarray arr = pyarray('d', [])

        # debug
        print self.__class__.__name__
        # end debug
        if self.baseptr_1.MatrixArray():
            for i in range(self.n_cols * self.n_rows):
                arr.append(self.baseptr_1.MatrixArray()[i])
            return arr
        else:
            raise ValueError("Can not get MatrixArray")

    @property
    def kind(self):
        return get_key(self.baseptr_1.Kind(), MatrixKindDict)


    @property
    def type(self):
        return get_key(self.baseptr_1.m2dType(), MatrixDict)

    #def void Add(self,size_t, void *):

    #def char * MatrixTypeString(self,MatrixType m):

    #def char * MatrixOutputString(self,MatrixType m):
