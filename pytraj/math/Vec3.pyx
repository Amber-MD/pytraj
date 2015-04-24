# distutils: language = c++
from cpython.array cimport array as pyarray
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr

from ..utils.check_and_assert import _import_numpy

# should we use numpy rather reinvent the wheel?
cdef class Vec3:
    def __cinit__(self, *args):
        cdef Vec3 vec
        cdef double x, y, z
        
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
        if self.thisptr is not NULL:
            del self.thisptr

    def __str__(self):
        x, y, z = self.tolist()
        return "Vec3 instance: %s %s %s" %(x, y, z)

    def Magnitude2(self):
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

    #def void operator /=(self,double xIn):
    def __idiv__(Vec3 self, double xIn):
        self.thisptr.divequal(xIn)
        return self

    #def Vec3 operator /(self,double xIn):
    def __div__(Vec3 self, double xIn):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr[0]/xIn
        return vec

    #def void operator *=(self,double xIn):
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
        # either return ndarray or memovryview
        return self.buffer1d[idx]

    def __setitem__(self, idx, value):
        if isinstance(value, (list, tuple)):
            value = pyarray('d', value)
        self.buffer1d[idx] = value

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
        has_np, np = _import_numpy()
        if has_np:
            np.asarray(self.buffer1d[:])

    @property
    def buffer1d(self):
        cdef double[:] arr = <double[:3]> self.thisptr.Dptr()
        return  arr
