# distutils: language = c++
from __future__ import absolute_import
from cython.view cimport array as cyarray
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from libcpp.vector cimport vector
from libc.math cimport sqrt
from libcpp.string cimport string
from cpython.array cimport array

from ..frame cimport Frame, _Frame
from ..core.box cimport _Box, Box

import numpy as np

__all__ = ['Grid', 'Matrix_3x3', 'distance_', 'torsion', ]

cdef class Grid:
    def __cinit__(self, *args):
        cdef size_t x, y, z

        self.thisptr = new _Grid[float]()

        if len(args) == 3:
            x, y, z = args
            self.resize(x, y, z)

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, idx):
        cdef size_t x, y, z
        x, y, z = idx
        return self.thisptr.element(x, y, z)

    def __setitem__(self, idx, value):
        cdef size_t x, y, z
        x, y, z = idx
        self.thisptr.setGrid(x, y, z, <float> value)

    @property
    def size(self):
        return self.thisptr.size()

    def resize(self, size_t x, size_t y, size_t z):
        self.thisptr.resize(x, y, z)

    @property
    def nx(self):
        return self.thisptr.NX()

    @property
    def ny(self):
        return self.thisptr.NY()

    @property
    def nz(self):
        return self.thisptr.NZ()

    def _element(self, int x, int y, int z):
        return self.thisptr.element(x, y, z)

    property data:
        def __get__(self):
            """return a copy of 3D array of Grid"""
            cdef size_t nx, ny, nz
            nx, ny, nz = self.nx, self.ny, self.nz
            cdef int i, j, k
            cdef cyarray carr = cyarray(shape=(nx, ny, nz),
                                        itemsize=sizeof(float), format="f")
            cdef float[:, :, :] myview = carr

            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        myview[i, j, k] = self.thisptr.element(i, j, k)
            return myview

    def to_ndarray(self):
        return np.asarray(self.data[:], dtype=np.float32)

    def tolist(self):
        return [[list(x) for x in y] for y in self.data]


cdef class Matrix_3x3:
    def __cinit__(self, Xin=None):
        """TODO:
             add doc
             Add: mat1 = Matrix_3x3(mat2)
                  (Cython complains "TypeError: 'src.Matrix_3x3.Matrix_3x3'
                   does not have the buffer interface")
        """
        cdef double[:] X
        cdef double _xin = 0.
        import numbers
        import itertools

        if Xin is None:
            self.thisptr = new _Matrix_3x3(0.0)
        elif isinstance(Xin, numbers.Number):
            _xin = Xin
            self.thisptr = new _Matrix_3x3(_xin)
        else:
            if len(Xin) == 9:
                # 1D
                X = np.asarray(Xin, dtype='f8')
            elif len(Xin) == 3:
                # 2D: 3x3
                _flat_list = []
                if hasattr(Xin, 'tolist'):
                    Xin = Xin.tolist()
                for x in Xin:
                    for x0 in x:
                        _flat_list.append(x0)
                X = np.asarray(_flat_list, dtype='f8')
            if X.shape[0] == 9:
                # Takes array of 9, row-major
                self.thisptr = new _Matrix_3x3(&X[0])
            elif X.shape[0] == 1:
                # Set all elements to the same number
                self.thisptr = new _Matrix_3x3(X[0])
            elif X.shape[0] == 3:
                # Set Set diagonal
                x, y, z = X
                self.thisptr = new _Matrix_3x3(x, y, z)
            else:
                raise ValueError("Must be array with length of None, 1, 3 or 9")

    def __iter__(self):
        cdef double[:] myview
        for myview in self.buffer2d[:]:
            yield myview

    def __getitem__(self, idx):
        return self.buffer2d[idx]

    def __setitem__(self, idx, value):
        self.buffer2d[idx] = value

    def __str__(self):
        txt = "Matrix 3x3: \n"
        for arr in self.buffer2d[:]:
            txt += " ".join(str(x) for x in arr)
            txt += "\n"
        return txt

    def __repr__(self):
        return self.__str__()

    def copy(self):
        cdef Matrix_3x3 other = Matrix_3x3()
        other.thisptr = new _Matrix_3x3(self.thisptr[0])
        return other

    def __dealloc__(self):
        """Free memory"""
        if self.thisptr:
            del self.thisptr

    def __imul__(Matrix_3x3 self, Matrix_3x3 other):
        """mat *= other"""

        self.thisptr[0].star_equal(other.thisptr[0])
        return self

    @property
    def row1(self):
        """
        Parameters: None
        ---------------

        Return
        ------
        Instance of Vec3
        """
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Row1()
        return vec

    @property
    def row2(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Row2()
        return vec

    @property
    def row3(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Row3()
        return vec

    @property
    def col1(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Col1()
        return vec

    @property
    def col2(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Col2()
        return vec

    @property
    def col3(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Col3()
        return vec

    def zeros(self):
        self.thisptr.Zero()

    def diagonalize(self, Vec3 vect):
        self.thisptr.Diagonalize(vect.thisptr[0])

    def diagonalize_sort(self, Vec3 vectds):
        self.thisptr.Diagonalize_Sort(vectds.thisptr[0])

    def diagonalize_sort_chirality(self, Vec3 vectds, int idx):
        self.thisptr.Diagonalize_Sort_Chirality(vectds.thisptr[0], idx)

    def transpose(self):
        self.thisptr.Transpose()

    def rotation_around_zaxis(self, idx, idy):
        self.thisptr.RotationAroundZ(idx, idy)

    def rotation_around_yaxis(self, idx, idz):
        self.thisptr.RotationAroundY(idx, idz)

    def __mul__(Matrix_3x3 self, other):
        cdef Matrix_3x3 mat_other, mat_result
        cdef Vec3 vec_other, vec_result

        if isinstance(other, Vec3):
            vec_other = other
            vec_result = Vec3()
            vec_result.thisptr[0] = self.thisptr[0] * vec_other.thisptr[0]
            return vec_result
        elif isinstance(other, Matrix_3x3):
            mat_other = other
            mat_result = Matrix_3x3()
            mat_result.thisptr[0] = self.thisptr[0] * mat_other.thisptr[0]
            return mat_result
        else:
            raise ValueError("Must be either Matrix_3x3 or Vec3")

    def calc_rotation_matrix(self, *args):
        """
        """
        cdef Vec3 vec
        cdef double theta
        cdef double x, y, z

        if len(args)== 2:
            vec, theta = args
            if not isinstance(vec, Vec3):
                raise ValueError("Must be a vector")
            self.thisptr.CalcRotationMatrix(vec.thisptr[0], theta)
        elif len(args) == 3:
            x, y, z = args
            self.thisptr.CalcRotationMatrix(x, y, z)
        else:
            raise ValueError('must be "Vec3, theta" or "x, y, z"')

    def rotation_angle(self):
        return self.thisptr.RotationAngle()

    def axis_of_rotation(self, theta):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.AxisOfRotation(theta)
        return vec

    def transpose_mult(self, Matrix_3x3 other):
        cdef Matrix_3x3 result
        result = Matrix_3x3()
        result.thisptr[0] = self.thisptr.TransposeMult(other.thisptr[0])
        return result

    def tolist(self):
        return [list(x) for x in self.buffer2d[:]]

    def to_ndarray(self, copy=True):
        if copy:
            return np.array(self.buffer2d)
        else:
            return np.asarray(self.buffer2d)

    def as_ndmatrix(self):
        return self.to_ndmatrix()

    def to_ndmatrix(self):
        """convert to numpy matrix as a memory view. No data copy is made"""
        mat = np.asmatrix(self.buffer2d[:])
        return mat

    property buffer2d:
        def __get__(self):
            cdef double[:, :] arr0 = <double[:3, :3]> self.thisptr.Dptr()
            return arr0

    property buffer1d:
        def __get__(self):
            cdef double[:] arr0 = <double[:9]> self.thisptr.Dptr()
            return arr0

    def __array__(self):
        return np.asarray(self.buffer2d)

cdef class Vec3:
    def __cinit__(self, *args):
        cdef Vec3 vec
        cdef double x, y, z

        self._own_memory = True
        if not args:
            self.thisptr = new _Vec3()
        else:
            if len(args) == 1 and isinstance(args[0], Vec3):
                vec = args[0]
                self.thisptr = new _Vec3(vec.thisptr[0])
            else:
                self.thisptr = new _Vec3()
                if isinstance(args[0], (list, tuple)):
                    x, y, z = args[0]
                    self.set_vec(x, y, z)
                else:
                    self.set_vec(*args)

    def __dealloc__(self):
        if self.thisptr is not NULL and self._own_memory:
            del self.thisptr

    def __str__(self):
        x, y, z = self.tolist()
        return "<Vec3: %s %s %s>" %(x, y, z)

    def __repr__(self):
        return self.__str__()

    def copy(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr = new _Vec3(self.thisptr[0])
        return vec

    def magnitude2(self):
        return self.thisptr.Magnitude2()

    def zeros(self):
        self.thisptr.Zero()

    def set_vec(self, *args):
        """args = either array of 3 elements or x, y, z"""
        cdef double x, y, z
        cdef double[:] X

        if len(args) == 1:
            # array
            X = args[0]
            if len(X) != 3:
                raise ValueError("must a an array of 3 elements")
            else:
                x = X[0]
                y = X[1]
                z = X[2]
                self.thisptr.SetVec(x, y, z)
        elif len(args) == 3:
            # x, y, z
            x, y, z = args
            self.thisptr.SetVec(x, y, z)
        else:
            raise ValueError("must a an array or x, y, z format")

    def normalize(self):
        return self.thisptr.Normalize()

    def angle(self, Vec3 othervec):
        return self.thisptr.Angle(othervec.thisptr[0])

    def assign(self, double[:] XYZ):
        self.thisptr.Assign(&XYZ[0])

    # def void operator /=(self,double xIn):
    def __idiv__(Vec3 self, double xIn):
        self.thisptr.divequal(xIn)
        return self

    # def Vec3 operator /(self,double xIn):
    def __div__(Vec3 self, double xIn):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr[0]/xIn
        return vec

    # def void operator *=(self,double xIn):
    # return "void": really?
    # Got: Segmentation fault (core dumped)
    def __imul__(Vec3 self, double xIn):
        self.thisptr.mulequal(xIn)
        return self

    def __mul__(Vec3 self, arg):
        cdef Vec3 rhs, vec
        cdef double xIn
        if isinstance(arg, Vec3):
            rhs = arg
            # return "double"
            return self.thisptr[0] * rhs.thisptr[0]
        else:
            # assuming "arg" is either double or int
            xIn = <double> arg
            vec = Vec3()
            vec.thisptr[0] = self.thisptr[0] * xIn
            return vec

    def __add__(Vec3 self, arg):
        cdef double xIn
        cdef Vec3 rhs, vec
        vec = Vec3()

        if isinstance(arg, Vec3):
            rhs = arg
            vec.thisptr[0]= self.thisptr[0] + rhs.thisptr[0]
            return vec
        else:
            xIn = arg
            vec.thisptr[0] = self.thisptr[0] + xIn
        return vec

    def __iadd__(Vec3 self, arg):
        cdef double xIn
        cdef Vec3 rhs
        if isinstance(arg, Vec3):
            rhs = arg
            self.thisptr.addequal(rhs.thisptr[0])
        else:
            xIn = arg
            self.thisptr.addequal(xIn)
        return self

    def __sub__(Vec3 self, arg):
        cdef double xIn
        cdef Vec3 rhs, vec
        vec = Vec3()

        if isinstance(arg, Vec3):
            rhs = arg
            vec.thisptr[0]= self.thisptr[0] - rhs.thisptr[0]
            return vec
        else:
            xIn = arg
            vec = self + (-xIn)
        return vec

    def __isub__(Vec3 self, arg):
        cdef double xIn
        cdef Vec3 rhs

        if isinstance(arg, Vec3):
            rhs = arg
            self.thisptr.subequal(rhs.thisptr[0])
        else:
            xIn = arg
            self.thisptr.addequal(xIn)
        return self

    def cross(self, Vec3 rhs):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Cross(rhs.thisptr[0])
        return vec

    def __getitem__(self, idx):
        return self.values[idx]

    def __setitem__(self, idx, value):
        self.values[idx] = value

    def is_zeros(self):
        return self.thisptr.IsZero()

    def neg(self):
        self.thisptr.Neg()

    def signed_angle(self, Vec3 v1, Vec3 v2):
        return self.thisptr.SignedAngle(v1.thisptr[0], v2.thisptr[0])

    def tolist(self):
        return list(self.buffer1d[:])

    def to_ndarray(self):
        """return a ndarray view of Vec3"""
        return np.asarray(self.buffer1d[:])

    @property
    def values(self):
        """return numpy array as memoryview of Vec3"""
        return self.to_ndarray()

    property buffer1d:
        def __get__(self):
            cdef double[:] arr = <double[:3]> self.thisptr.Dptr()
            return arr

cdef extern from "DistRoutines.h" nogil:
    ctypedef enum ImagingType:
        NOIMAGE=0
        ORTHO
        NONORTHO

    double DIST2_ImageNonOrtho "DIST2_ImageNonOrtho"(const _Vec3 &, const _Vec3 &, const _Matrix_3x3 &, const _Matrix_3x3 &)
    double DIST2_ImageNonOrthoRecip(const _Vec3 &, const _Vec3 &, double, int*, const _Matrix_3x3 &)
    double DIST2_ImageOrtho(const _Vec3 &, const _Vec3 &, const _Box &)
    double DIST2_NoImage_from_ptr "DIST2_NoImage"(const double*, const double*)
    double DIST2_NoImage(const _Vec3 &, const _Vec3 &)
    double DIST_NoImage "DIST2_NoImage"(const _Vec3 &, const _Vec3 &)
    double DIST2(const double*, const double*, ImagingType, const _Box &,
                 const _Matrix_3x3 &, const _Matrix_3x3 &)


def distance_(double[:, :, :] p):
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
# distutils: language = c++


cdef class ImagedAction:
    def __cinit__(self):
        self.thisptr = new _ImagedAction()

    def __dealloc__(self):
        del self.thisptr

    # def ImagedAction(self):

    # def void InitImaging(self,bint imageIn):

    # def void SetupImaging(self,Box::BoxType parmboxtype):

    # def bint ImagingEnabled(self):

    # def bint UseImage(self):

    # def ImagingType ImageType(self):

# distutil: language = c++
import math

cdef extern from "TorsionRoutines.h" nogil:
    # create alias to avoid: ambiguous overloaded method
    double C_Torsion "Torsion" (const double *, const double *, const double *, const double *)
    double Pucker_AS(const double*, const double*, const double*, const double*,
                     const double*, double &)
    double Pucker_CP(const double*, const double*, const double*, const double*,
                     const double*, const double*, int, double &, double &)
    double C_CalcAngle "CalcAngle" (const double*, const double*, const double*)


def torsion(double[:, :, :] p):
    cdef double[:] out = np.empty(p.shape[0])
    cdef int i

    if p.shape[1] != 4 or p.shape[2] != 3:
        raise ValueError("shape of input array must be (n_frames, 4, 3)")

    for i in range(p.shape[0]):
        out[i] = math.degrees(C_Torsion(&p[i, 0, 0], &p[i, 1, 0],
                                         &p[i, 2, 0], &p[i, 3, 0]))
    return np.asarray(out)


def angle(double[:, :, :] p):
    cdef double[:] out = np.empty(p.shape[0])
    cdef int i

    if p.shape[1] != 3 or p.shape[2] != 3:
        raise ValueError("shape of input array must be (n_frames, 4, 3)")

    for i in range(p.shape[0]):
        out[i] = math.degrees(C_CalcAngle(&p[i, 0, 0], &p[i, 1, 0],
                                           &p[i, 2, 0]))
    return np.asarray(out)
