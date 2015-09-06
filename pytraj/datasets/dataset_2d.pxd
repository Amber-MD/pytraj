# distutils: language = c++
from __future__ import absolute_import
from libcpp.vector cimport vector
from .base cimport _DataSet, DataSet, DataType
from ..math.Matrix cimport *
from .dataset_2d cimport _DataSet_2D, DataSet_2D


cdef extern from "DataSet_2D.h": 
    # DataSet_2D.h
    ctypedef enum MatrixType "DataSet_2D::MatrixType":
        pass
    ctypedef enum MatrixKind "DataSet_2D::MatrixKind":
        pass
    cdef cppclass _DataSet_2D "DataSet_2D" (_DataSet):
        _DataSet_2D() 
        _DataSet_2D(DataType tIn, int wIn, int pIn)
        # virtual methods
        int Allocate2D(size_t, size_t) 
        int AllocateHalf(size_t) 
        int AllocateTriangle(size_t) 
        double GetElement(size_t, size_t) const  
        size_t Nrows() const  
        size_t Ncols() const  
        double * MatrixArray() const  
        MatrixKind Kind "MatrixKind"() const  
        # end virtual methods

        void Add(size_t, const void *)
        const char * MatrixTypeString(MatrixType m)
        const char * MatrixOutputString(MatrixType m)

cdef class DataSet_2D (DataSet):
    cdef _DataSet_2D* baseptr_1


#ctypedef Matrix[double].iterator iterator
ctypedef vector[double] Darray

cdef extern from "DataSet_MatrixDbl.h": 
    cdef cppclass _DatasetMatrixDouble "DataSet_MatrixDbl" (_DataSet_2D):
        _DatasetMatrixDouble() 
        double& index_opr "operator[]"(size_t idx)
        @staticmethod
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate2D(size_t x, size_t y)
        int AllocateHalf(size_t x)
        int AllocateTriangle(size_t x)
        #void Write2D(_CpptrajFile&, int, int) const 
        double GetElement(size_t x, size_t y) const 
        size_t Nrows() const 
        size_t Ncols() const 
        #double * MatrixArray() const # not implemented
        MatrixKind Kind() const 
        # make alias to avoid naming conflict with DataSet (DataType)
        MatrixType matType "Type"() const 
        unsigned int Nsnapshots() const 
        void IncrementSnapshots() 
        double& Element(size_t x, size_t y)
        int AddElement(double d)
        void SetElement(size_t x, size_t y, double d)
        #iterator begin() 
        #iterator end() 
        const Darray& Vect() const 
        Darray& V1() 
        void AllocateVector(size_t vsize)
        #Darray.iterator v1begin() 
        #Darray.iterator v1end() 
        void SetTypeAndKind(MatrixType tIn, MatrixKind kIn)
        void StoreMass(const Darray& mIn)
        const Darray& Mass() const 


cdef class DatasetMatrixDouble (DataSet_2D):
    cdef _DatasetMatrixDouble* thisptr
    cdef bint py_free_mem


cdef extern from "DataSet_MatrixFlt.h": 
    cdef cppclass _DatasetMatrixFloat  "DataSet_MatrixFlt" (_DataSet_2D):
        _DataSet_MatrixFlt() 
        float& index_opr "operator[]" (size_t idx)
        @staticmethod
        _DataSet * Alloc() 


cdef class DatasetMatrixFloat(DataSet_2D):
    cdef _DatasetMatrixFloat * thisptr
    cdef bint py_free_mem
