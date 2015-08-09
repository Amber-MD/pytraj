# distutils: language = c++


cdef class ReplicaDimArray:
    def __cinit__(self):
        self.thisptr = new _ReplicaDimArray()

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, int idx):
        return self.thisptr.index_opr(idx)

    def AddRemdDimension(self,int d):
        self.thisptr.AddRemdDimension(d)

    def clear(self):
        self.thisptr.clear()

    def Ndims(self):
        return self.thisptr.Ndims()

    def Description(self,int idx):
        return self.thisptr.Description(idx)

    #def bint operator !=(self, ReplicaDimArray rhs):
    def __richcmp__(self, ReplicaDimArray rsh, opt):
        pass


