# distutils: language = c++


cdef class Traj_AmberRestart:
    def __cinit__(self):
        self.thisptr = new _Traj_AmberRestart()

    def __dealloc__(self):
        del self.thisptr

    def Traj_AmberRestart(self):

    def BaseIOtype * Alloc(self):

