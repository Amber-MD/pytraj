# distutil: language = c++

from ..core.Box cimport _Box, Box
from pytraj.math.Vec3 cimport _Vec3, Vec3
from pytraj.math.Matrix_3x3 cimport _Matrix_3x3, Matrix_3x3
from libc.math cimport sqrt
from libcpp.string cimport string
from cpython.array cimport array
from ..Frame cimport Frame, _Frame

cdef extern from "DistRoutines.h" nogil:
    ctypedef enum ImagingType:
        NOIMAGE=0 
        ORTHO
        NONORTHO
    
    double DIST2_ImageNonOrtho "DIST2_ImageNonOrtho"(const _Vec3 &, const _Vec3 &, const _Matrix_3x3 &, const _Matrix_3x3 &)
    double DIST2_ImageNonOrthoRecip(const _Vec3 &, const _Vec3&, double, int*, const _Matrix_3x3&)
    double DIST2_ImageOrtho(const _Vec3&, const _Vec3&, const _Box &)
    double DIST2_NoImage_from_ptr "DIST2_NoImage"(const double*, const double*)
    double DIST2_NoImage( const _Vec3&, const _Vec3& )
    double DIST_NoImage "DIST2_NoImage"( const _Vec3&, const _Vec3& )
    double DIST2(const double*, const double*, ImagingType, const _Box &, 
                 const _Matrix_3x3&, const _Matrix_3x3&)


def distance_(double[:, :, :] p):
   import numpy as np
   cdef double[:] out = np.empty(p.shape[0])
   cdef int i

   if p.shape[1] != 2 or p.shape[2] != 3:
       raise ValueError("shape of input array must be (n_frames, 2, 3)")

   for i in range(p.shape[0]):
       out[i] = DIST2_NoImage_from_ptr(&p[i, 0, 0], &p[i, 1, 0])
   return np.sqrt(np.asarray(out))

def dist2_image_nonOrtho(Vec3 v1, Vec3 v2, Matrix_3x3 m1, Matrix_3x3 m2):
    return DIST2_ImageNonOrtho(v1.thisptr[0], v2.thisptr[0], m1.thisptr[0], m2.thisptr[0])

def dist_noimage(Vec3 v1, Vec3 v2):
    return DIST_NoImage(v1.thisptr[0], v2.thisptr[0])

def distance(p1, p2, image=None, image_type=None, *args, **kwd):
    cdef Vec3 v1
    cdef Vec3 v2
    cdef double[:] arr1
    cdef double[:] arr2
    cdef int i

    if not image:
        if isinstance(p1, Vec3) and isinstance(p2, Vec3):
            v1 = <Vec3> p1
            v2 = <Vec3> p2
        else:
            v1 = Vec3(p1[0], p1[1], p1[2])
            v2 = Vec3(p2[0], p2[1], p2[2])
        return sqrt(DIST2_NoImage(v1.thisptr[0], v2.thisptr[0]))
    else:
        raise NotImplementedError("not yet supported")
