# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector
#from pytraj.MaskToken cimport *
from pytraj.cpp_vector cimport vector as cppvector


cdef extern from "AtomMask.h": 
    cdef cppclass _AtomMask "AtomMask":
        _AtomMask()
        _AtomMask(const string&)
        _AtomMask(int, int)
        _AtomMask(int)
        _AtomMask(vector[int], int)
        _AtomMask(const _AtomMask &)
        #_AtomMask & operator =(const _AtomMask &)
        const vector [int]& Selected()const 
        cppvector[int].const_iterator begin()const 
        cppvector[int].const_iterator end()const 
        int back()const 
        int Nselected()const 
        const int & index_opr "operator[]"(int idx)const 
        const char * MaskString()const 
        const string& MaskExpression()const 
        bint MaskStringSet()const 
        bint None()const 
        bint IsCharMask()const 
        void ResetMask()
        void ClearSelected()
        void InvertMask() except +
        int NumAtomsInCommon(const _AtomMask&)
        void AddSelectedAtom(int i)
        void AddAtom(int)
        void AddAtoms(const vector [int]&)
        void AddAtomRange(int, int)
        void AddMaskAtPosition(const _AtomMask&, int)
        void PrintMaskAtoms(const char *)const 
        #int SetMaskString(const char *)
        int SetMaskString(const string&)
        void SetupIntMask(const char *, int, int)
        void SetupCharMask(const char *, int, int)
        bint AtomInCharMask(int)const 
        bint AtomsInCharMask(int, int)const 
        void SetNatom(int a)
        int ConvertToCharMask()
        int ConvertToIntMask()
        void MaskInfo()const 
        void BriefMaskInfo()const 
        #inline token_iterator begintoken()const 
        #inline token_iterator endtoken()const 

#ctypedef fused charstring:
#    char*
#    string

cdef class AtomMask:
    cdef _AtomMask* thisptr
