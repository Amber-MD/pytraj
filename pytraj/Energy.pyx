# distutils: language = c++


cdef class Energy_Amber:
    def __cinit__(self):
        self.thisptr = new _Energy_Amber()

    def __dealloc__(self):
        del self.thisptr

    def E_bond(self, Frame frame, Topology top, AtomMask mask):
        return self.thisptr.E_bond(frame.thisptr[0], top.thisptr[0], mask.thisptr[0])

    def E_angle(self, Frame frame, Topology top, AtomMask mask):
        return self.thisptr.E_angle(frame.thisptr[0], top.thisptr[0], mask.thisptr[0])

    def E_torsion(self, Frame frame, Topology top, AtomMask mask):
        return self.thisptr.E_torsion(frame.thisptr[0], top.thisptr[0], mask.thisptr[0])

    def E_bond(self, Frame frame, Topology top, AtomMask mask):
        return self.thisptr.E_bond(frame.thisptr[0], top.thisptr[0], mask.thisptr[0])

    def E_14_nonbond(self, Frame frame, Topology top, AtomMask mask, double Eq14):
        return self.thisptr.E_14_Nonbond(frame.thisptr[0], top.thisptr[0], mask.thisptr[0], Eq14)

    def E_nonbond(self, Frame frame, Topology top, AtomMask mask, double Eelec):
        return self.thisptr.E_Nonbond(frame.thisptr[0], top.thisptr[0], mask.thisptr[0], Eelec)

    def set_debug(self, int d):
        self.thisptr.SetDebug(d)

    def print_timing(self):
        self.thisptr.PrintTiming()
