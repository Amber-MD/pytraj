# distutil: language = c++

from pytraj.core.Box cimport _Box, Box
from pytraj.math.Vec3 cimport _Vec3, Vec3
from pytraj.math.Matrix_3x3 cimport _Matrix_3x3, Matrix_3x3
from libc.math cimport sqrt

cdef extern from "DistRoutines.h":
    ctypedef enum ImagingType:
        NOIMAGE=0 
        ORTHO
        NONORTHO
    
    double DIST2_ImageNonOrtho "DIST2_ImageNonOrtho"(const _Vec3 &, const _Vec3 &, const _Matrix_3x3 &, const _Matrix_3x3 &)
    double DIST2_ImageNonOrthoRecip(const _Vec3 &, const _Vec3&, double, int*, const _Matrix_3x3&)
    double DIST2_ImageOrtho(const _Vec3&, const _Vec3&, const _Box &)
    double DIST2_NoImage(const double*, const double*)
    double DIST2_NoImage( const _Vec3&, const _Vec3& )
    double DIST_NoImage "DIST2_NoImage"( const _Vec3&, const _Vec3& )
    double DIST2(const double*, const double*, ImagingType, const _Box &, 
                 const _Matrix_3x3&, const _Matrix_3x3&)

