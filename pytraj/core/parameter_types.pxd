# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string

ctypedef vector[_BondParmType] BondParmArray
ctypedef vector[_BondType] BondArray
ctypedef vector[_AngleParmType] AngleParmArray
ctypedef vector[_AngleType] AngleArray
ctypedef vector[_DihedralType] DihedralArray
ctypedef vector[_HB_ParmType] HB_ParmArray
ctypedef vector[_NonbondType] NonbondArray
ctypedef vector[_LES_AtomType] LES_Array
ctypedef vector[_CmapGridType] CmapGridArray
ctypedef vector[_CmapType] CmapArray

# should I need to define fused type?
ctypedef fused ptype:
    _BondParmType
    _BondType
    _AngleParmType

cdef extern from "ParameterTypes.h": 
    cdef cppclass _AngleType "AngleType":
        _AngleType() 
        _AngleType(int a1, int a2, int a3, int idx)
        inline int A1() const 
        inline int A2() const 
        inline int A3() const 
        inline int Idx() const 

    cdef cppclass _NonbondParmType "NonbondParmType":
        _NonbondParmType() 
        _NonbondParmType(int n, const vector[int]& nbi, const NonbondArray& nba, const HB_ParmArray& hba)
        inline bint HasNonbond() const 
        inline int Ntypes() const 
        const vector[int]& NBindex() const 
        const NonbondArray& NBarray() const 
        const HB_ParmArray& HBarray() const 
        const _NonbondType& NBarray(int i) const 
        const _HB_ParmType& HBarray(int i) const 
        int GetLJindex(int type1, int type2) const 


    cdef cppclass _LES_AtomType "LES_AtomType":
        _LES_AtomType() 
        _LES_AtomType(int t, int c, int i)
        inline int Type() const 
        inline int Copy() const 
        inline int ID() const 


    cdef cppclass _AngleParmType "AngleParmType":
        _AngleParmType() 
        _AngleParmType(double tk, double teq)
        inline double Tk() const 
        inline double Teq() const 


    cdef cppclass _CmapType "CmapType":
        _CmapType() 
        _CmapType(int a1, int a2, int a3, int a4, int a5, int i)
        inline int A1() const 
        inline int A2() const 
        inline int A3() const 
        inline int A4() const 
        inline int A5() const 
        inline int Idx() const 


    cdef cppclass _LES_ParmType "LES_ParmType":
        _LES_ParmType() 
        _LES_ParmType(int na, int nt, const vector[double]& fac)
        inline bint HasLES() const 
        inline int Ntypes() const 
        inline int Ncopies() const 
        const vector[double]& FAC() const 
        const LES_Array& Array() const 
        void SetTypes(int n, const vector[double]& f)
        void AddLES_Atom(const _LES_AtomType& lat)


    cdef cppclass _HB_ParmType "HB_ParmType":
        _HB_ParmType() 
        _HB_ParmType(double a, double b, double c)
        inline double Asol() const 
        inline double Bsol() const 
        inline double HBcut() const 


    cdef cppclass _NonbondType "NonbondType":
        _NonbondType() 
        _NonbondType(double a, double b)
        inline double A() const 
        inline double B() const 


    cdef cppclass _BondParmType "BondParmType":
        _BondParmType() 
        _BondParmType(double rk, double req)
        inline double Rk() const 
        inline double Req() const 


    cdef cppclass _CmapGridType "CmapGridType":
        _CmapGridType() 
        _CmapGridType(int r, const vector[double]& g)
        inline int Resolution() const 
        inline const vector[double]& Grid() const 


    # ParameterTypes.h
    ctypedef enum Dtype "DihedralType::Dtype":
        NORMAL "DihedralType::NORMAL"
        IMPROPER "DihedralType::IMPROPER"
        END "DihedralType::END"
        BOTH "DihedralType::BOTH"

    cdef cppclass _DihedralType "DihedralType":
        _DihedralType() 
        _DihedralType(int a1, int a2, int a3, int a4, int idx)
        _DihedralType(int a1, int a2, int a3, int a4, Dtype t, int i)
        int A1()
        int A2()
        int A3()
        int A4()
        Dtype Type()
        int Idx()


    cdef cppclass _BondType "BondType":
        _BondType() 
        _BondType(int a1, int a2, int idx)
        inline int A1() const 
        inline int A2() const 
        inline int Idx() const 
        void SetIdx(int i)


    cdef cppclass _CapParmType "CapParmType":
        _CapParmType() 
        _CapParmType(int n, double c, double x, double y, double z)
        inline bint HasWaterCap() const 
        inline int NatCap() const 
        inline double CutCap() const 
        inline double xCap() const 
        inline double yCap() const 
        inline double zCap() const 


cdef class AngleType:
    cdef _AngleType* thisptr

cdef class NonbondParmType:
    cdef _NonbondParmType* thisptr

cdef class LES_AtomType:
    cdef _LES_AtomType* thisptr

cdef class AngleParmType:
    cdef _AngleParmType* thisptr

cdef class CmapType:
    cdef _CmapType* thisptr

cdef class LES_ParmType:
    cdef _LES_ParmType* thisptr

cdef class HB_ParmType:
    cdef _HB_ParmType* thisptr

cdef class NonbondType:
    cdef _NonbondType* thisptr

cdef class BondParmType:
    cdef _BondParmType* thisptr

cdef class CmapGridType:
    cdef _CmapGridType* thisptr

cdef class DihedralType:
    cdef _DihedralType* thisptr

cdef class BondType:
    cdef _BondType* thisptr

cdef class CapParmType:
    cdef _CapParmType* thisptr
