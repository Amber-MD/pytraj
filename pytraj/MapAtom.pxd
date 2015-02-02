# distutils: language = c++
from pytraj.Atom cimport *


cdef extern from "MapAtom.h": 
    cdef cppclass _MapAtom "MapAtom" (_Atom):
        _MapAtom() 
        _MapAtom(const _MapAtom&)
        _MapAtom(const _Atom&)
        #_MapAtom& operator =(const _MapAtom&)
        bint IsChiral() const 
        bint BoundToChiral() const 
        bint IsMapped() const 
        bint Complete() const 
        bint IsUnique() const 
        const string& _AtomID() const 
        const string& Unique() const 
        int Nduplicated() const 
        char CharName() const 
        void IsDuplicated() 
        void SetMapped() 
        void SetComplete() 
        void SetChiral() 
        void SetBoundToChiral() 
        void SetAtomID(const string&)
        void SetUnique(const string&)
        void SetNotMapped() 
        void SetNotComplete() 
        void SetNotChiral() 

cdef class MapAtom (Atom):
    # _Atom* thisptr
    cdef _MapAtom* thisptr_ma
