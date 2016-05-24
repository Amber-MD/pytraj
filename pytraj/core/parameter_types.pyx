# distutils: language = c++

from cpython.array cimport array as pyarray

cdef class AngleType:
    def __cinit__(self, arg=None):
        cdef int a1, a2, a3, idx
        if not arg:
            self.thisptr = new _AngleType()
        else:
            if len(arg) == 4:
                a1, a2, a3, idx = arg
                self.thisptr = new _AngleType(a1, a2, a3, idx)
            else:
                raise ValueError()

    def __dealloc__(self):
        del self.thisptr

    @property
    def idx(self):
        return self.thisptr.Idx()

    @property
    def indices(self):
        """return atom indices as a python array"""
        # cdef pyarray arr = pyarray('i', [self.thisptr.Idx(), self.thisptr.A1(),
        #                                 self.thisptr.A2(), self.thisptr.A3()])
        cdef pyarray arr = pyarray('i', [self.thisptr.A1(),
                                         self.thisptr.A2(), self.thisptr.A3()])
        return arr


cdef class LES_AtomType:
    def __cinit__(self):
        self.thisptr = new _LES_AtomType()

    def __dealloc__(self):
        del self.thisptr

    # def LES_AtomType(self):

    # def LES_AtomType(self,int t, int c, int i):

    @property
    def type(self):
        return self.thisptr.Type()

    def copy(self):
        return self.thisptr.Copy()

    @property
    def ID(self):
        return self.thisptr.ID()

cdef class AngleParmType:
    def __cinit__(self, arg=None):
        cdef double tk, teq
        if not arg:
            self.thisptr = new _AngleParmType()
        else:
            if len(arg) == 2:
                tk, teq = arg
                self.thisptr = new _AngleParmType(tk, teq)

    def __dealloc__(self):
        del self.thisptr

    @property
    def Tk(self):
        return self.thisptr.Tk()

    @property
    def Teq(self):
        return self.thisptr.Teq()

cdef class BondParmType:
    def __cinit__(self, arg=None):
        cdef double rk, req
        if not arg:
            self.thisptr = new _BondParmType()
        else:
            if len(arg) == 2:
                rk, req = arg
                self.thisptr = new _BondParmType(rk, req)
            else:
                raise ValueError()

    def __dealloc__(self):
        del self.thisptr

    @property
    def Rk(self):
        return self.thisptr.Rk()

    @property
    def Req(self):
        return self.thisptr.Req()

cdef class DihedralType:
    def __cinit__(self, arg=None):
        cdef int a1, a2, a3, a4, idx
        cdef Dtype dhtype
        if not arg:
            self.thisptr = new _DihedralType()
        else:
            # if len(arg) == 5:
            #    a1, a2, a3, a4, idx = arg
            #    self.thisptr = new _DihedralType(a1, a2, a3, a4, idx)
            # elif len(arg) == 6:
            #    a1, a2, a3, a4, idx, dhtype = arg
            #    self.thisptr = new _DihedralType(a1, a2, a3, a4, dhtype, idx)
            # else:
            #    raise ValueError()
            raise ValueError()

    def __dealloc__(self):
        del self.thisptr

    @property
    def type(self):
        return self.thisptr.Type()

    @property
    def idx(self):
        return self.thisptr.Idx()

    @property
    def indices(self):
        """return atom indices as a python array"""
        # cdef pyarray arr = pyarray('i', [self.thisptr.Idx(), self.thisptr.A1(), self.thisptr.A2(),
        #                                 self.thisptr.A3(), self.thisptr.A4()])
        cdef pyarray arr = pyarray('i', [self.thisptr.A1(), self.thisptr.A2(),
                                         self.thisptr.A3(), self.thisptr.A4()])
        return arr

cdef class BondType:
    def __cinit__(self, arg=None):
        cdef int a1, a2, idx
        if not arg:
            self.thisptr = new _BondType()
        else:
            if len(arg) == 3:
                a1, a2, idx = arg
                self.thisptr = new _BondType(a1, a2, idx)
            else:
                raise ValueError()

    def __dealloc__(self):
        del self.thisptr

    property idx:
        def __get__(self):
            return self.thisptr.Idx()

        def __set__(self, int i):
            self.thisptr.SetIdx(i)

    @property
    def indices(self):
        """return atom indices as a python array"""
        # cdef pyarray arr = pyarray('i', [self.thisptr.Idx(), self.thisptr.A1(),
        # self.thisptr.A2()])
        cdef pyarray arr = pyarray('i', [self.thisptr.A1(), self.thisptr.A2()])
        return arr


cdef class DihedralParmType:
    def __cinit__(self, arg=None):
        cdef double k, n, p, e, b
        if not arg:
            self.thisptr = new _DihedralParmType()
        else:
            if len(arg) == 2:
                k, p = arg
                self.thisptr = new _DihedralParmType(k, p)
            elif len(arg) == 5:
                k, n, p, e, b = arg
                self.thisptr = new _DihedralParmType(k, n, p, e, b)
            else:
                raise ValueError()

    def __dealloc__(self):
        del self.thisptr

    @property
    def Pk(self):
        return self.thisptr.Pk()

    @property
    def Pn(self):
        return self.thisptr.Pn()

    @property
    def phase(self):
        return self.thisptr.Phase()

    property SCEE:
        def __get__(self):
            return self.thisptr.SCEE()

        def __set__(self, double s):
            self.thisptr.SetSCEE(s)

    property SCNB:
        def __get__(self):
            return self.thisptr.SCNB()

        def __set__(self, double s):
            self.thisptr.SetSCNB(s)
