# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from pytraj.cpptraj_dict import BoxTypeDict, get_key


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
        txt = "<Box instance with x, y, z, alpha, beta, gamma = %s>" % boxlisttxt
        return txt

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, idx):
        """add fancy indexing?"""
        #return self.thisptr.index_opr(idx)
        return self.tolist()[idx]

    def __setitem__(self, idx, double value):
        cdef double* ptr
        if not isinstance(idx, (long, int)):
            raise NotImplementedError("support only integer indexing, not slice")
        ptr = &(self.thisptr.index_opr(idx))
        ptr[0] = value

    def __iter__(self):
        for i in range(6):
            yield self[i]

    @classmethod
    def get_all_boxtypes(cls):
        return [x.lower() for x in BoxTypeDict.keys()]

    @classmethod
    def help(cls):
        print cls.get_all_boxtypes()


    @property
    def bname(self):
        return self.thisptr.TypeName().decode()
    
    def set_beta_lengths(self, double beta, double xin, double yin, double zin):
        self.thisptr.SetBetaLengths(beta, xin, yin, zin)

    def set_box(self, double[:] boxIn):
        self.thisptr.SetBox(&boxIn[0])

    def set_trunc_oct(self):
        self.thisptr.SetTruncOct()

    def set_nobox(self):
        self.thisptr.SetNoBox()

    def set_missing_info(self, Box boxinst):
        self.thisptr.SetMissingInfo(boxinst.thisptr[0])

    def to_recip(self,Matrix_3x3 ucell, Matrix_3x3 recip):
        return self.thisptr.ToRecip(ucell.thisptr[0], recip.thisptr[0])

    @property
    def btype(self):
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

    def tolist(self):
        cdef int i
        cdef vector[double] v
        cdef double* ptr = self.thisptr.boxPtr()

        for i in range(6):
            v.push_back(deref(ptr+i))
        return v
