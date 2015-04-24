# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.MapAtom cimport MapAtom, _MapAtom
from pytraj.Topology cimport _Topology, Topology


cdef extern from "AtomMap.h": 
    cdef cppclass _AtomMap "AtomMap":
        _AtomMap() 
        _MapAtom& operator[](int)
        const _MapAtom& operator[](int i) const 
        int Natom() const 
        void SetDebug(int d)
        int Setup(const _Topology&)
        int SetupResidue(const _Topology&, int)
        void ResetMapping() 
        bint BondIsRepeated(int, int) const 
        void Determine_AtomIDs() 
        void MarkAtomComplete(int, bint)
        void CheckForComplete_Atoms() 
        #int SymmetricAtoms(const _Topology&, _AtomIndexArray&, int)


cdef class AtomMap:
    cdef _AtomMap* thisptr
