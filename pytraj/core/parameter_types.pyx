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

cdef class NonbondParmType:

    def __cinit__(self, arg=None):
        """TODO: add and check"""
        cdef int n
        cdef vector[int] nbi
        cdef NonbondArray nba
        cdef HB_ParmArray hba

        if not arg:
            self.thisptr = new _NonbondParmType()
        else:
            if len(arg) == 4:
                n, nbi, nba_, hba_ = arg
            # to be added
            #nba = convert_objlist_to_vector(BondType, nba_, NonbondArray)
            #hba = convert_objlist_to_vector(HB_ParmType, hba_, HB_ParmArray)
            # self.thisptr = new _NonbondParmType(n, nbi, nba, hba)

    def __dealloc__(self):
        del self.thisptr

    def has_nonbond(self):
        return self.thisptr.HasNonbond()

    @property
    def n_types(self):
        return self.thisptr.Ntypes()

    @property
    def NBindex(self):
        return self.thisptr.NBindex()

    # def  NonbondArray NBarray(self):
    #    """return list of NonbondType instances?"""
    #    cdef NonbondArray na
    #    cdef NonbondType nbt
    #    cdef _NonbondType nbt_
    #    na = self.thisptr.NBarray()
    #    nbt_list = []
    #    for nbt_ in na:
    #        nbt = NonbondType()
    #        del nbt.thisptr
    #        # to be added

    # def  HB_ParmArray HBarray(self):

    # def  NonbondType NBarray(self,int i):

    # def  HB_ParmType HBarray(self,int i):

    def GetLJindex(self, int type1, int type2):
        return self.thisptr.GetLJindex(type1, type2)

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

cdef class CmapType:
    def __cinit__(self, arg):
        cdef int a1, a2, a3, a4, a5, i
        if not arg:
            self.thisptr = new _CmapType()
        else:
            if len(arg) == 6:
                a1, a2, a3, a4, a5, i = arg
                self.thisptr = new _CmapType(a1, a2, a3, a4, a5, i)

    def __dealloc__(self):
        del self.thisptr

    @property
    def idx(self):
        return self.thisptr.Idx()

    @property
    def indices(self):
        """return atom indices as a python array"""
        # cdef pyarray arr = pyarray('i', [self.thisptr.Idx(), self.thisptr.A1(),
        #                                 self.thisptr.A2(), self.thisptr.A3(),
        #                                 self.thisptr.A4(), self.thisptr.A5()])
        cdef pyarray arr = pyarray('i', [self.thisptr.A1(),
                                         self.thisptr.A2(), self.thisptr.A3(),
                                         self.thisptr.A4(), self.thisptr.A5()])
        return arr


cdef class HB_ParmType:
    def __cinit__(self, arg=None):
        cdef double a, b, c
        if not arg:
            self.thisptr = new _HB_ParmType()
        else:
            if len(arg) == 3:
                a, b, c = arg
                self.thisptr = new _HB_ParmType(a, b, c)

    def __dealloc__(self):
        del self.thisptr

    @property
    def Asol(self):
        return self.thisptr.Asol()

    @property
    def Bsol(self):
        return self.thisptr.Bsol()

    @property
    def HBcut(self):
        return self.thisptr.HBcut()

cdef class NonbondType:
    def __cinit__(self, arg=None):
        if not arg:
            self.thisptr = new _NonbondType()
        else:
            if len(arg) == 2:
                a, b = arg
                self.thisptr = new _NonbondType(a, b)
            else:
                raise ValueError()

    def __dealloc__(self):
        del self.thisptr

    @property
    def A(self):
        return self.thisptr.A()

    @property
    def B(self):
        return self.thisptr.B()

cdef class ChamberParmType:
    def __cinit__(self):
        self.thisptr = new _ChamberParmType()

    def __dealloc__(self):
        del self.thisptr

    def HasChamber(self):
        return self.thisptr.HasChamber()

    def has_cmap(self):
        return self.thisptr.HasCmap()

    @property
    def FF_version(self):
        return self.thisptr.FF_Version()

    @property
    def FF_type(self):
        return self.thisptr.FF_Type()

    # def  BondArray UB(self):

    # def  BondParmArray UBparm(self):

    # def  DihedralArray Impropers(self):

    # def  DihedralParmArray ImproperParm(self):

    # def  NonbondArray LJ14(self):

    # def  CmapGridArray CmapGrid(self):

    # def  CmapArray Cmap(self):

    # def SetLJ14(self, NonbondArray nb):
    #    self.thisptr.SetLJ14(nb)

    def SetChamber(self, int i, string s):
        self.thisptr.SetChamber(i, s)

    # def SetUB(self, BondArray ub, BondParmArray ubp):
    #    self.thisptr.SetUB(ub, ubp)

    def SetImproper(self, im_py, imp_py):
        cdef DihedralArray im
        cdef DihedralParmArray imp
        cdef DihedralType im_inst
        cdef DihedralParmType imp_inst

        # convert list of im_py to vector of Dihedral
        for im_inst in im_py:
            im.push_back(im_inst.thisptr[0])
        for imp_inst in imp_py:
            imp.push_back(imp_inst.thisptr[0])
        self.thisptr.SetImproper(im, imp)

    def AddCmapGrid(self, CmapGridType g):
        self.thisptr.AddCmapGrid(g.thisptr[0])

    def AddCmapTerm(self, CmapType c):
        self.thisptr.AddCmapTerm(c.thisptr[0])

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

cdef class CmapGridType:
    def __cinit__(self, *args):
        cdef int r
        cdef vector[double] g
        if not args:
            self.thisptr = new _CmapGridType()
        else:
            if len(args) == 2:
                r, g = args
                self.thisptr = new _CmapGridType(r, g)

    def __dealloc__(self):
        del self.thisptr

    @property
    def resolution(self):
        return self.thisptr.Resolution()

    @property
    def grid(self):
        return self.thisptr.Grid()

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


cdef class CapParmType:
    def __cinit__(self, arg=None):
        cdef int n
        cdef double c, x, y, z
        if not arg:
            self.thisptr = new _CapParmType()
        else:
            if len(arg) == 5:
                n, c, x, y, z = arg
                self.thisptr = new _CapParmType(n, c, x, y, z)
            else:
                raise ValueError()

    def __dealloc__(self):
        del self.thisptr

    def has_water_cap(self):
        return self.thisptr.HasWaterCap()

    @property
    def NatCap(self):
        return self.thisptr.NatCap()

    @property
    def cutcap(self):
        return self.thisptr.CutCap()

    @property
    def x_cap(self):
        return self.thisptr.xCap()

    @property
    def y_cap(self):
        return self.thisptr.yCap()

    @property
    def z_cap(self):
        return self.thisptr.zCap()

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
