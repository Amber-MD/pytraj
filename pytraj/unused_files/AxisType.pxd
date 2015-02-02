# distutils: language = c++
from pycpptraj.Topology cimport *


cdef extern from "AxisType.h": 
    # AxisType.h
    ctypedef enum NAType "NA_Base::NAType":
        UNKNOWN_BASE "NA_Base::UNKNOWN_BASE"
        ADE "NA_Base::ADE"
        CYT "NA_Base::CYT"
        GUA "NA_Base::GUA"
        THY "NA_Base::THY"
        URA "NA_Base::URA"
    # AxisType.h
    ctypedef enum PSType "NA_Base::PSType":
        PHOS "NA_Base::PHOS"
        O4p "NA_Base::O4p"
        C1p "NA_Base::C1p"
        C2p "NA_Base::C2p"
        C3p "NA_Base::C3p"
        C4p "NA_Base::C4p"
    cdef cppclass _NA_Base "NA_Base":
        _NA_Base() 
        _NA_Base(const _NA_Base&)
        #_NA_Base& operator =(const _NA_Base&)
        NAType ID_BaseFromName(const _NameType&)
        _NA_Base(const _Topology&, int, NAType)
        void SetInput_Frame(const _Frame&)
        void Print_AtomNames() const 
        NAType Type() const 
        int ResNum() const 
        char BaseChar() const 
        const _Frame& Ref() const 
        const _Frame& Input() const 
        const _AtomMask& InputFitMask() const 
        const _AtomMask& RefFitMask() const 
        const char * _AtomName(int i) const 
        bint HasPatom() const 
        bint HasO4atom() const 
        bint HasSugar_Atoms() const 
        const char * ResName() const 
        const char * RefName(int i) const 
        int HBidx(int i) const 
        const double * HBxyz(int i) const 
        const double * Pxyz() const 
        const double * O4xyz() const 
        const double * C1xyz() const 
        const double * C2xyz() const 
        const double * C3xyz() const 
        const double * C4xyz() const 
    cdef cppclass _NA_Axis "NA_Axis":
        _NA_Axis() 
        void SetupBaseAxis(const _Matrix_3x3&, const _Vec3&, int)
        _NA_Axis(int, int, bint)
        void StoreRotMatrix(const _Matrix_3x3&, const _Vec3&)
        void PrintAxisInfo(const char *) const 
        void FlipYZ() 
        void FlipXY() 
        const _Matrix_3x3& Rot() const 
        const _Vec3& Oxyz() const 
        const _Vec3& Rx() const 
        const _Vec3& Ry() const 
        const _Vec3& Rz() const 
        int Res1() const 
        int Res2() const 
        bint IsAnti() const 
