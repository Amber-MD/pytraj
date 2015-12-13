# distutils: language = c++
from ..math.cpp_math cimport _Matrix_3x3, Matrix_3x3, _Vec3, Vec3
from ..cpp_vector cimport vector


cdef extern from "Box.h": 
    ctypedef enum BoxType "Box::BoxType":
        NOBOX "Box::NOBOX"
        ORTHO "Box::ORTHO"
        TRUNCOCT "Box::TRUNCOCT"
        RHOMBIC "Box::RHOMBIC"
        NONORTHO "Box::NONORTHO"
    cdef cppclass _Box "Box":
        _Box() 
        _Box(const double *)
        _Box(const _Box &)
        _Box(_Matrix_3x3 const)
        #_Box & operator =(const _Box &)
        const char * TypeName() const 
        void SetBetaLengths(double, double, double, double)
        void SetBox(const double *)
        void SetTruncOct() 
        void SetNoBox() 
        void SetMissingInfo(const _Box &)
        double ToRecip(_Matrix_3x3 &, _Matrix_3x3 &)const 
        void SetX(double xin)
        void SetY(double yin)
        void SetZ(double zin)
        void SetAlpha(double ain)
        void SetBeta(double bin)
        void SetGamma(double gin)
        BoxType Type() const 
        double BoxX() const 
        double BoxY() const 
        double BoxZ() const 
        double Alpha() const 
        double Beta() const 
        double Gamma() const 
        bint HasBox() const 
        _Vec3 Center() const 
        _Vec3 Lengths() const 
        double * boxPtr() 
        #const double * boxPtr() const 
        #const double& index_opr "operator[]"(int idx)const 
        double& index_opr "operator[]"(int idx)

cdef class Box:
    cdef _Box* thisptr
