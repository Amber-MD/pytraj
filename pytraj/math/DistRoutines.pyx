# distutil: language = c++
from libcpp.string cimport string
from cpython.array cimport array
from pytraj.Frame cimport Frame, _Frame

#from pytraj.Vec3 cimport _Vec3, Vec3

def dist2_image_nonOrtho(Vec3 v1, Vec3 v2, Matrix_3x3 m1, Matrix_3x3 m2):
    return DIST2_ImageNonOrtho(v1.thisptr[0], v2.thisptr[0], m1.thisptr[0], m2.thisptr[0])

def dist_noimage(Vec3 v1, Vec3 v2):
    return DIST_NoImage(v1.thisptr[0], v2.thisptr[0])

#DIST2_ImageNonOrthoRecip(const _Vec3 &, const _Vec3&, double, int*, const _Matrix_3x3&)
#DIST2_ImageOrtho(const _Vec3&, const _Vec3&, const _Box &)
#DIST2_NoImage(const double*, const double*)
#DIST2_NoImage( const _Vec3&, const _Vec3& )
#DIST_NoImage( const _Vec3&, const _Vec3& )
#DIST2(const double*, const double*, ImagingType, const _Box &, 
#      const _Matrix_3x3&, const _Matrix_3x3&)

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

def distance_frames(Frame f1, Frame f2, image=False, image_type="", *args, **kwd):
    # TODO : addd *args and **kwd
    if f1.n_atoms != f2.n_atoms:
        raise ValueError("two frames must have the same number of atoms")
    image_type = image_type.encode()
    return array('d', _distance_frames(f1.thisptr[0], f2.thisptr[0], image, image_type, f1.n_atoms))

cdef double[:] _distance_frames(_Frame f1, _Frame f2, bint image, 
                                string image_type, int natoms):
    # TODO : 
    #     + extend this method
    #     + test cases
    cdef int i
    cdef double[:] arr0 = array('d', [-1]*natoms)

    if not image:
        for i in range(natoms):
            arr0[i] = sqrt(DIST2_NoImage(f1.XYZ(i), f2.XYZ(i)))
    else:
        raise NotImplementedError("not yet supported for image")
    return arr0
