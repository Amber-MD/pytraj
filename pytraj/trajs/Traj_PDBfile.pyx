# distutils: language = c++


cdef class Traj_PDBfile:
    def __cinit__(self):
        self.thisptr = new _Traj_PDBfile()

    def __dealloc__(self):
        del self.thisptr

    def alloc(self):
        cdef BaseIOtype baseio = BaseIOtype()
        baseio.baseptr0 = self.thisptr.Alloc()
        return baseio

    #@classmethod
    #def WriteHelp(cls):
    #    _Traj_PDBfile.WriteHelp()

