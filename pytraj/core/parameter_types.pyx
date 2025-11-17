# distutils: language = c++

import numpy as np

cdef class AngleType:
    def __cinit__(self):
        self.thisptr = new _AngleType()

    def __dealloc__(self):
        del self.thisptr

    @property
    def idx(self):
        return self.thisptr.Idx()

    @property
    def indices(self):
        """return atom indices as a python array"""
        return np.array([self.thisptr.A1(),
                         self.thisptr.A2(), self.thisptr.A3()])

cdef class BondType:
    def __cinit__(self):
        self.thisptr = new _BondType()

    def __dealloc__(self):
        del self.thisptr

    @property
    def idx(self):
        return self.thisptr.Idx()

    @property
    def indices(self):
        """return atom indices as a python array"""
        return np.array([self.thisptr.A1(), self.thisptr.A2()])

cdef class DihedralType:
    def __cinit__(self):
        self.thisptr = new _DihedralType()

    def __dealloc__(self):
        del self.thisptr

    @property
    def idx(self):
        return self.thisptr.Idx()

    @property
    def indices(self):
        """return atom indices as a python array"""
        return np.array([self.thisptr.A1(), self.thisptr.A2(),
                         self.thisptr.A3(), self.thisptr.A4()])