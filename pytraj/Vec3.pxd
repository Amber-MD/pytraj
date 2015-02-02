# distutils: language = c++
from libcpp.vector cimport vector


cdef extern from "Vec3.h": 
    cdef cppclass _Vec3 "Vec3":
        _Vec3() 
        _Vec3(const _Vec3& rhs)
        _Vec3(double vx, double vy, double vz)
        _Vec3(double vxyz)
        _Vec3(const double * XYZ)
        _Vec3(const float * XYZ)
        _Vec3(const int * XYZ)
        #_Vec3& operator =(const _Vec3& rhs)
        void Assign(const double * XYZ)
        void divequal "operator /=" (double xIn)
        _Vec3 operator /(double xIn) const 
        void mulequal "operator *=" (double xIn)
        _Vec3 operator *(double xIn) const 
        void addequal "operator +=" (double xIn)
        _Vec3 operator +(double xIn) const 
        void subequal "operator -=" (const _Vec3& rhs)
        _Vec3 operator -(const _Vec3& rhs) const 
        void addequal "operator +=" (const _Vec3& rhs)
        _Vec3 operator +(const _Vec3& rhs) const 
        double operator *(const _Vec3& rhs) const 
        _Vec3 Cross(const _Vec3& rhs) const 
        #double operator[](int idx) const 
        #double& index_opr "operator[]"(int idx)
        double& index_opr "operator[]"(int idx) const
        double Magnitude2() const 
        void Zero() 
        bint IsZero() const 
        void Neg() 
        void SetVec(double vx, double vy, double vz)
        double Normalize() 
        void Print(const char *) const 
        double Angle(const _Vec3&) const 
        double SignedAngle(const _Vec3&, const _Vec3&) const 
        #const double * Dptr() const 
        double * Dptr() 


cdef class Vec3:
    cdef _Vec3* thisptr

