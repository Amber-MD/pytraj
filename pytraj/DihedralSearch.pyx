# distutils: language = c++


cdef class DihedralSearch:
    def __cinit__(self):
        self.thisptr = new _DihedralSearch()

    def __dealloc__(self):
        del self.thisptr

    #def mask_it begin(self):
    #def mask_it end(self):

    def n_dihedrals(self):
        return self.thisptr.Ndihedrals()

    def list_known_types(self):
        self.thisptr.ListKnownTypes()

    def offset_help(self):
        self.thisptr.OffsetHelp()

    def get_type(self, string typein):
        return self.thisptr.GetType(typein)

    def search_for(self, DihedralType2 dhtype):
        return self.thisptr.SearchFor(dhtype)

    def search_for_args(self, ArgList arglist):
        self.thisptr.SearchForArgs(arglist.thisptr[0])

    def search_for_new_type(self, int off, string an0, string an1, string an2, string an3, string name):
        return self.SearchForNewType(off, an0, an1, an2, an3, name)

    def search_for_all(self):
        return self.thisptr.SearchForAll()

#    def find_dihedrals(self, Topology top, Range r):
#        return self.thisptr.FindDihedrals(top.thisptr[0], r.thisptr[0])

    def clear(self):
        self.thisptr.Clear()

    def clear_found(self):
        self.thisptr.ClearFound()

    def print_types(self):
        self.thisptr.PrintTypes()

    def moving_atoms(self, Topology top, int atom0, int atom1):
        cdef AtomMask atm = AtomMask()
        atm.thisptr[0] = self.thisptr.MovingAtoms(top.thisptr[0],atom0, atom1)
        return atm

# private class of DihedralSearch
# Dont' use
#cdef class DihedralMask:
#    def __cinit__(self, *args):
#        """DihedralMask() 
#        or DihedralMask(a0, a1, a2, a3, res, name, dhtype)
#        """
#        cdef int a0, a1, a2, a3, res
#        cdef string name
#        cdef DihedralType dhtype
#        
#        if not args:
#            self.thisptr = new _DihedralMask()
#        elif len(args) == 7:
#            a0, a1, a2, a3, res, name, dhtype = args
#            self.thisptr = new _DihedralMask(a0, a1, a2, a3, res, name, dhtype)
#
#    def __dealloc__(self):
#        del self.thisptr
#
#    def a0(self):
#        return self.thisptr.A0()
#
#    def a1(self):
#        return self.thisptr.A1()
#
#    def a2(self):
#        return self.thisptr.A2()
#
#    def a3(self):
#        return self.thisptr.A3()
#
#    def res_num(self):
#        return self.thisptr.ResNum()
#
#    def name(self):
#        return self.thisptr.Name()
#
#    def none(self):
#        return self.thisptr.None()
#
#    def type(self):
#        return self.thisptr.Type()
#
#
