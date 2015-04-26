# distutils: language = c++
from __future__ import absolute_import
from cython cimport view
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarray
from pytraj.cpptraj_dict import BoxTypeDict, get_key
from pytraj.utils import _import_numpy


cdef class Box:
    def __cinit__(self, *args):
        cdef double[:] boxIn 
        cdef Box rhs

        if not args:
            self.thisptr = new _Box()
        elif len(args) == 1:
            if isinstance(args[0], Box):
                rhs = args[0]
                self.thisptr = new _Box(rhs.thisptr[0])
            else:
                boxIn = args[0]
                self.thisptr = new _Box(&boxIn[0])
        else: 
            raise ValueError("")

    def __str__(self):
        boxlisttxt = ", ".join([str(tmp) for tmp in self.tolist()])
        txt = "<Box with x, y, z, alpha, beta, gamma = %s>" % boxlisttxt
        return txt

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, idx):
        """add fancy indexing?"""
        return self.data[idx]

    def __setitem__(self, idx, value):
        self.data[idx] = value

    def __iter__(self):
        for x in self.data:
            yield x

    @classmethod
    def all_box_types(cls):
        return [x.lower() for x in BoxTypeDict.keys()]

    @property
    def name(self):
        return self.thisptr.TypeName().decode()
    
    def set_beta_lengths(self, double beta, double xin, double yin, double zin):
        self.thisptr.SetBetaLengths(beta, xin, yin, zin)

    def set_box_from_array(self, boxIn):
        # try to cast array-like to python array
        # list, tuple are ok too
        cdef pyarray arr0 = pyarray('d', boxIn)
        cdef double[:] myview = arr0
        
        self.thisptr.SetBox(&myview[0])

    def set_trunc_oct(self):
        self.thisptr.SetTruncOct()

    def set_nobox(self):
        self.thisptr.SetNoBox()

    def set_missing_info(self, Box boxinst):
        self.thisptr.SetMissingInfo(boxinst.thisptr[0])

    def to_recip(self,Matrix_3x3 ucell, Matrix_3x3 recip):
        return self.thisptr.ToRecip(ucell.thisptr[0], recip.thisptr[0])

    @property
    def type(self):
        return get_key(self.thisptr.Type(), BoxTypeDict).lower()

    property x:
        def __get__(self):
            return self.thisptr.BoxX()
        def __set__(self, double value):
            self.thisptr.SetX(value)

    property y:
        def __get__(self):
            return self.thisptr.BoxY()
        def __set__(self, double value):
            self.thisptr.SetY(value)

    property z:
        def __get__(self):
            return self.thisptr.BoxZ()
        def __set__(self, double value):
            self.thisptr.SetZ(value)

    property alpha:
        def __get__(self):
            return self.thisptr.Alpha()
        def __set__(self, double value):
            self.thisptr.SetAlpha(value)

    property beta:
        def __get__(self):
            return self.thisptr.Beta()
        def __set__(self, double value):
            self.thisptr.SetBeta(value)

    property gamma:
        def __get__(self):
            return self.thisptr.Gamma()
        def __set__(self, double value):
            self.thisptr.SetGamma(value)

    def has_box(self):
        return self.thisptr.HasBox()

    @property
    def center(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Center()
        return vec

    @property
    def lengths(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Lengths()
        return vec

    @property
    def data(self):
        """memoryview for box array"""
        cdef double[:] arr0
        arr0 = <double[:6]> self.thisptr.boxPtr()
        return arr0

    def tolist(self):
        return list(self.data[:])

    def to_ndarray(self):
        has_np, np = _import_numpy()
        if not has_np:
            raise ImportError("need numpy")
        else:
            return np.asarray(self.data[:])
