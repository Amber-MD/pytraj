# distutils: language = c++


cdef class DataSet_GridFlt:
    def __cinit__(self):
        self.baseptr0 = <_DataSet*> new _DataSet_GridFlt()
        self.baseptr_1 = <_DataSet_3D*> self.baseptr0
        self.thisptr = <_DataSet_GridFlt*> self.baseptr0
        self.py_free_mem = True

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    #def float operator[](self,size_t idx):
    def __getitem__(self, size_t idx):
        return self.thisptr.index_opr(idx)

    def Alloc(self):
        # cpptraj: return DataSet*
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.Alloc()
        return dset

    #def Grid[float] InternalGrid(self):

    @property
    def size(self):
        return self.thisptr.Size()

    ##def int Sync(self):

    def Info(self):
        self.thisptr.Info()

    def Allocate3D(self,size_t x, size_t y, size_t z):
        return self.thisptr.Allocate3D(x, y, z)

    #def void Write3D(self,CpptrajFile, int, int, int):

    def GetElement(self,int x, int y, int z):
        return self.thisptr.GetElement(x, y, z)

    def SetElement(self,int x, int y, int z, float v):
        self.thisptr.SetElement(x, y, z, v)

    def NX(self):
        return self.thisptr.NX()

    def NY(self):
        return self.thisptr.NY()

    def NZ(self):
        return self.thisptr.NZ()

    #def iterator begin(self):

    #def iterator end(self):

    #def  long int Increment(self, Vec3, float):

    #def  long int Increment(self, double *, float):

    #def  long int Increment(self,int, int, int, float):

    #def float GridVal(self,int x, int y, int z):

    #def long int CalcIndex(self,int i, int j, int k):

