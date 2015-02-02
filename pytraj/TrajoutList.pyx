# distutils: language = c++


cdef class TrajoutList:
    def __cinit__(self):
        self.thisptr = new _TrajoutList()

    def __dealloc__(self):
        del self.thisptr

    #def void Clear(self):

    #def void SetDebug(self,int):

    #def int AddEnsembleTrajout(self, ArgList, TopologyList, int):

    #def int AddTrajout(self, ArgList, TopologyList):

    #def int WriteTrajout(self,int, Topology *, Frame *):

    #def void CloseTrajout(self):

    #def void List(self):

    #def bint Empty(self):

    #def ArgIt argbegin(self):

    #def ArgIt argend(self):

