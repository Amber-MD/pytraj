# distutils: language = c++


cdef class DataSet_Vector:
    def __cinit__(self):
        self.thisptr = new _DataSet_Vector()

    def __dealloc__(self):
        del self.thisptr

    #cdef DataSet * Alloc(self):

    @property
    def Size(self):
        return self.thisptr.Size()

    def Sync(self):
        return self.thisptr.Sync()

    def Info(self):
        self.thisptr.Info()

    #def int Allocate1D(self,size_t):

    #def  void Add(self,size_t, void *):

    #def double Dval(self,size_t):

    #def double Xcrd(self,size_t idx):

    #def void WriteBuffer(self,CpptrajFile, size_t):

    #def void SetIred(self):

    #def bint IsIred(self):

    #def void reset(self):

    #def void Resize(self,size_t s):

    #def void Resize(self,size_t s, Vec3 v):

    #def bint Empty(self):

    #def  Vec3 operator[](self,int i):

    #def Vec3 operator[](self,int i):

    #def  Vec3 OXYZ(self,int i):

    #def void ReserveVecs(self,size_t n):

    #def void AddVxyz(self, Vec3 v):

    #def void AddVxyz(self, Vec3 v, Vec3 c):

    ##def const_iterator begin(self):

    ##def const_iterator end(self):

    #def Vec3 Back(self):

    #def int CalcSphericalHarmonics(self,int):

    #def ComplexArray SphericalHarmonics(self,int):

    #def double SphericalHarmonicsNorm(self,int):

