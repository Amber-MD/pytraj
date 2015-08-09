# distutils: language = c++
from pytraj.core.NameType cimport _NameType, NameType

cdef extern from "Residue.h": 
    cdef cppclass _Residue "Residue":
        _Residue()
        _Residue(int onum, const _NameType& resname, int first_AtomIn)
        inline void SetLastAtom(int i)
        inline void SetOriginalNum(int i)
        inline int FirstAtom() const 
        inline int LastAtom() const 
        inline int OriginalResNum() const 
        inline const char * c_str() const 
        inline const _NameType& Name() const 
        inline int NumAtoms() const 
        inline bint NameIsSolvent() const 

cdef class Residue:
    cdef _Residue* thisptr

