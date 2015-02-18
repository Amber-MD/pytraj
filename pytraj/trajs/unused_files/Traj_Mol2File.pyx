# distutils: language = c++


cdef class Traj_Mol2File:
    def __cinit__(self):
        self.thisptr = new _Traj_Mol2File()

    def __dealloc__(self):
        del self.thisptr

    #def Traj_Mol2File(self):

    #def BaseIOtype * Alloc(self):

