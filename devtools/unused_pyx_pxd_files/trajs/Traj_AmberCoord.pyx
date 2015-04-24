# distutils: language = c++


cdef class Traj_AmberCoord (TrajectoryIO):
    def __cinit__(self):
        self.baseptr_1 = <_TrajectoryIO*> new _Traj_AmberCoord()
        self.thisptr = <_Traj_AmberCoord*> self.baseptr_1

    def __dealloc__(self):
        del self.thisptr

    def alloc(self):
        cdef BaseIOtype baseio = BaseIOtype()
        baseio.baseptr0 = self.thisptr.Alloc()
        return baseio

    # turn off: why do we need this?
    #def write_help(self):
    #    self.thisptr.WriteHelp()
