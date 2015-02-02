# distutils: language = c++
from libcpp.vector cimport vector


cdef extern from "ReplicaDimArray.h": 
    # ReplicaDimArray.h
    ctypedef enum RemDimType "ReplicaDimArray::RemDimType":
        UNKNOWN "ReplicaDimArray::UNKNOWN"
        TEMPERATURE "ReplicaDimArray::TEMPERATURE"
        PARTIAL "ReplicaDimArray::PARTIAL"
        HAMILTONIAN "ReplicaDimArray::HAMILTONIAN"
        PH "ReplicaDimArray::PH"
    cdef cppclass _ReplicaDimArray "ReplicaDimArray":
        _ReplicaDimArray() 
        int index_opr "operator[]"(int idx) const 
        void AddRemdDimension(int d)
        void clear() 
        int Ndims() const 
        const char * Description(int idx) const 
        bint operator ! =(const _ReplicaDimArray& rhs) const 


cdef class ReplicaDimArray:
    cdef _ReplicaDimArray* thisptr

