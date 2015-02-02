# distutils: language = c++
from libcpp.string cimport string
from pytraj.NameType cimport *


cdef extern from "MaskToken.h": 
    # MaskToken.h
    ctypedef enum MaskTokenType "MaskToken::MaskTokenType":
        OP_NONE "MaskToken::OP_NONE"
        ResNum "MaskToken::ResNum"
        ResName "MaskToken::ResName"
        AtomNum "MaskToken::AtomNum"
        AtomName "MaskToken::AtomName"
        AtomType "MaskToken::AtomType"
        AtomElement "MaskToken::AtomElement"
        SelectAll "MaskToken::SelectAll"
        OP_AND "MaskToken::OP_AND"
        OP_OR "MaskToken::OP_OR"
        OP_NEG "MaskToken::OP_NEG"
        OP_DIST "MaskToken::OP_DIST"
    cdef cppclass _MaskToken "MaskToken":
        _MaskToken() 
        const char* MaskTypeString[]
        const char * TypeName() const 
        void Print() const 
        int SetToken(MaskTokenType, const string&)
        int SetDistance(string&)
        void SetOperator(MaskTokenType)
        inline MaskTokenType Type() const 
        inline int Res1() const 
        inline int Res2() const 
        inline const _NameType& Name() const 
        inline bint OnStack() const 
        inline bint Within() const 
        inline bint ByAtom() const 
        inline double Distance() const 
        void SetOnStack() 


cdef class MaskToken:
    cdef _MaskToken* thisptr

