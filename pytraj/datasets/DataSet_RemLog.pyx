# distutils: language = c++


cdef class DataSet_RemLog:
    def __cinit__(self):
        self.thisptr = new _DataSet_RemLog()

    def __dealloc__(self):
        del self.thisptr

    #def DataSet_RemLog(self):

    #def DataSet * Alloc(self):

    #def void AllocateReplicas(self,int):

    #def void AddRepFrame(self,int rep, ReplicaFrame frm):

    #def  ReplicaFrame RepFrame(self,int exch, int rep):

    #def int NumExchange(self):

    #def bint ValidEnsemble(self):

    #def void TrimLastExchange(self):

    #def size_t Size(self):

    #def int Sync(self):

    #def void Info(self):

    #def void Add(self,size_t, void *):

cdef class ReplicaFrame:
    def __cinit__(self):
        self.thisptr = new _ReplicaFrame()

    def __dealloc__(self):
        del self.thisptr

    #def int SetTremdFrame(self, char *, TmapType):

    #def int SetHremdFrame(self, char *, vector[int]):

    #def int ReplicaIdx(self):

    #def int PartnerIdx(self):

    #def int CoordsIdx(self):

    #def bint Success(self):

    #def double Temp0(self):

    #def double PE_X1(self):

    #def double PE_X2(self):

