#void  distutils: language = c++


cdef extern from "Molecule.h": 
    cdef cppclass _Molecule "Molecule":
        _Molecule()
        _Molecule(int begin, int end)
        void SetFirst(int begin)
        void SetLast(int last)
        void SetSolvent() 
        void SetNoSolvent() 
        inline int BeginAtom() const 
        inline int EndAtom() const 
        inline bint IsSolvent() const 
        inline int NumAtoms() const 

cdef class Molecule:
    cdef _Molecule* thisptr

