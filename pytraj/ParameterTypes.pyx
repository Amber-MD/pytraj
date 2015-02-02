# distutils: language = c++

# TODO: 
#     + make convert_objlist_to_vector works right
#        from cy_utils cimport convert_objlist_to_vector
#     + uncomment some methods
# 

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
    def A1(self):
        return self.thisptr.A1()

    @property
    def A2(self):
        return self.thisptr.A2()

    @property
    def A3(self):
        return self.thisptr.A3()

    @property
    def Idx(self):
        return self.thisptr.Idx()

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
            #self.thisptr = new _NonbondParmType(n, nbi, nba, hba)

    def __dealloc__(self):
        del self.thisptr

    def HasNonbond(self):
        return self.thisptr.HasNonbond()

    @property
    def Ntypes(self):
        return self.thisptr.Ntypes()

    @property 
    def NBindex(self):
        return self.thisptr.NBindex()

    #def  NonbondArray NBarray(self):
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

    #def  HB_ParmArray HBarray(self):

    #def  NonbondType NBarray(self,int i):

    #def  HB_ParmType HBarray(self,int i):

    def GetLJindex(self,int type1, int type2):
        return self.thisptr.GetLJindex(type1, type2)

cdef class LES_AtomType:
    def __cinit__(self):
        self.thisptr = new _LES_AtomType()

    def __dealloc__(self):
        del self.thisptr

    #def LES_AtomType(self):

    #def LES_AtomType(self,int t, int c, int i):

    @property
    def Type(self):
        return self.thisptr.Type()

    def Copy(self):
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
    def A1(self):
        return self.thisptr.A1()

    @property
    def A2(self):
        return self.thisptr.A2()

    @property
    def A3(self):
        return self.thisptr.A3()

    @property
    def A4(self):
        return self.thisptr.A4()

    @property
    def A5(self):
        return self.thisptr.A5()

    @property
    def Idx(self):
        return self.thisptr.Idx()

# not need LES_ParmType now
#cdef class LES_ParmType:
#    def __cinit__(self):
#        self.thisptr = new _LES_ParmType()
#
#    def __dealloc__(self):
#        del self.thisptr
#
#    def LES_ParmType(self):
#
#    def LES_ParmType(self,int na, int nt, vector[double] fac):
#
#    def  bint HasLES(self):
#
#    def  int Ntypes(self):
#
#    def  int Ncopies(self):
#
#    def  vector[double] FAC(self):
#
#    def  LES_Array Array(self):
#
#    def void SetTypes(self,int n, vector[double] f):
#
#    def void AddLES_Atom(self, LES_AtomType lat):

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

    def HasCmap(self):
        return self.thisptr.HasCmap()

    @property
    def FF_Version(self):
        return self.thisptr.FF_Version()

    @property
    def FF_Type(self):
        return self.thisptr.FF_Type()

    #def  BondArray UB(self):

    #def  BondParmArray UBparm(self):

    #def  DihedralArray Impropers(self):

    #def  DihedralParmArray ImproperParm(self):

    #def  NonbondArray LJ14(self):

    #def  CmapGridArray CmapGrid(self):

    #def  CmapArray Cmap(self):

    #def SetLJ14(self, NonbondArray nb):
    #    self.thisptr.SetLJ14(nb)

    def SetChamber(self,int i, string s):
        self.thisptr.SetChamber(i, s)

    #def SetUB(self, BondArray ub, BondParmArray ubp):
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
        cdef int a1, a2, a3, a4, idx, i
        cdef Dtype t
        if not arg:
            self.thisptr = new _DihedralType()
        else:
            if len(arg) == 5:
                a1, a2, a3, a4, idx = arg
                self.thisptr = new _DihedralType(a1, a2, a3, a4, idx)
            elif len(arg) == 6:
                a1, a2, a3, a4, t, i = arg
                self.thisptr = new _DihedralType(a1, a2, a3, a4, t, i)
            else:
                raise ValueError()

    def __dealloc__(self):
        del self.thisptr

    @property
    def A1(self):
        return self.thisptr.A1()

    @property
    def A2(self):
        return self.thisptr.A2()

    @property
    def A3(self):
        return self.thisptr.A3()

    @property
    def A4(self):
        return self.thisptr.A4()

    @property
    def type(self):
        return self.thisptr.Type()

    @property
    def idx(self):
        return self.thisptr.Idx()

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

    @property
    def A1(self):
        return self.thisptr.A1()

    @property
    def A2(self):
        return self.thisptr.A2()

    @property
    def Idx(self):
        return self.thisptr.Idx()

    @Idx.setter
    def Idx(self,int i):
        self.thisptr.SetIdx(i)


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
    def CutCap(self):
        return self.thisptr.CutCap()

    @property
    def xCap(self):
        return self.thisptr.xCap()

    @property
    def yCap(self):
        return self.thisptr.yCap()

    @property
    def zCap(self):
        return self.thisptr.zCap()

cdef class DihedralParmType:
    def __cinit__(self, arg=None):
        cdef double k, n, p, e, b
        if not arg:
            self.thisptr = new _DihedralParmType()
        else:
            if len(arg) == 2:
                k ,p = arg
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
    def Phase(self):
        return self.thisptr.Phase()

    @property
    def SCEE(self):
        return self.thisptr.SCEE()

    @SCEE.setter
    def SCEE(self,double s):
        self.thisptr.SetSCEE(s)

    @property
    def SCNB(self):
        return self.thisptr.SCNB()

    @SCNB.setter
    def SCNB(self,double s):
        self.thisptr.SetSCNB(s)
