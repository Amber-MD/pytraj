# distutils: language = c++
# don't know why getting import error if not include "pytraj.math"
from pytraj.math.Vec3 cimport _Vec3, Vec3


cdef extern from "Matrix_3x3.h": 
    cdef cppclass _Matrix_3x3 "Matrix_3x3":
        _Matrix_3x3() 
        _Matrix_3x3(const _Matrix_3x3&)
        _Matrix_3x3(const double *)
        _Matrix_3x3(double)
        _Matrix_3x3(double, double, double)
        #_Matrix_3x3& operator =(const _Matrix_3x3&)
        double operator[](int idx) const 
        double& operator[](int idx)
        _Vec3 Row1() 
        _Vec3 Row2() 
        _Vec3 Row3() 
        _Vec3 Col1() 
        _Vec3 Col2() 
        _Vec3 Col3() 
        void Zero() 
        void Print(const char *) const 
        int Diagonalize(_Vec3&)
        int Diagonalize_Sort(_Vec3&)
        int Diagonalize_Sort_Chirality(_Vec3&, int)
        void Transpose() 
        _Matrix_3x3& star_equal "operator *=" (const _Matrix_3x3&)
        void RotationAroundZ(double, double)
        void RotationAroundY(double, double)
        void CalcRotationMatrix(const _Vec3&, double)
        void CalcRotationMatrix(double, double, double)
        double RotationAngle() 
        _Vec3 AxisOfRotation(double)
        _Vec3 operator *(const _Vec3& rhs) const 
        #_Vec3 TransposeMult(const _Vec3& rhs) const  #not yet implemented in cpptraj?
        _Matrix_3x3 operator *(const _Matrix_3x3&) const 
        _Matrix_3x3 TransposeMult(const _Matrix_3x3&) const 
        #const double * Dptr() const 
        double * Dptr() 

cdef class Matrix_3x3:
    cdef _Matrix_3x3* thisptr
