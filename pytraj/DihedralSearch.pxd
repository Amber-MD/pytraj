# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.Topology cimport *
from pytraj.datasets.DataSet cimport *
from pytraj.Range cimport *
from pytraj.ArgList cimport *


cdef extern from "DihedralSearch.h": 
    # Note: _DihedralMask is private class of DihedralSearch
    # Don't use
    cdef cppclass _DihedralMask "DihedralSearch::DihedralMask":
        _DihedralMask() 
        _DihedralMask(int, int, int, int, int, const string&, DihedralType2)
        int A0() const 
        int A1() const 
        int A2() const 
        int A3() const 
        int ResNum() const 
        const string& Name() const 
        bint None() const 
        DihedralType2 Type() const 

    # DihedralSearch.h
    ctypedef enum DihedralType2 "DihedralSearch::DihedralType":
        PHI "DihedralSearch::PHI"
        PSI "DihedralSearch::PSI"
        CHIP "DihedralSearch::CHIP"
        OMEGA "DihedralSearch::OMEGA"
        ALPHA "DihedralSearch::ALPHA"
        BETA "DihedralSearch::BETA"
        GAMMA "DihedralSearch::GAMMA"
        DELTA "DihedralSearch::DELTA"
        EPSILON "DihedralSearch::EPSILON"
        ZETA "DihedralSearch::ZETA"
        NU1 "DihedralSearch::NU1"
        NU2 "DihedralSearch::NU2"
        CHIN "DihedralSearch::CHIN"
        NDIHTYPE "DihedralSearch::NDIHTYPE"

    cdef cppclass _DihedralSearch "DihedralSearch":
        #mask_it begin() const 
        #mask_it end() const 
        int Ndihedrals() const 
        _DihedralSearch() 
        void ListKnownTypes() 
        void OffsetHelp() 
        DihedralType2 GetType(const string&)
        int SearchFor(DihedralType2)
        void SearchForArgs(_ArgList&)
        int SearchForNewType(int, const string&, const string&, const string&, const string&, const string&)
        int SearchForAll() 
        int FindDihedrals(const _Topology&, const _Range&)
        void Clear() 
        void ClearFound() 
        void PrintTypes() 
        _AtomMask MovingAtoms(const _Topology&, int, int)


    # private class?
    #cdef cppclass DihedralToken "DihedralSearch::DihedralToken":
    #    DihedralToken() 
    #    DihedralToken(int, const _NameType&, const _NameType&, const _NameType&, const _NameType&, const string&)
    #    DihedralToken(const DIH_TYPE&, DihedralType2)
    #    DihedralMask FindDihedral_Atoms(const _Topology&, int)
    #    const string& Name() const 
    #    DihedralType2 Type() const 
    #    void Set_AtomName(int i, const _NameType& n)


cdef class DihedralSearch:
    cdef _DihedralSearch* thisptr

#cdef class DihedralMask: 
#    cdef _DihedralMask* thisptr
