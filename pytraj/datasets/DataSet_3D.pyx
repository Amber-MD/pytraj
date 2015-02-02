# distutils: language = c++


cdef class DataSet_3D (DataSet):
    def __cinit__(self):
        self.baseptr_1 = <_DataSet_3D*> self.baseptr0

    def __dealloc__(self):
        # since this is ABC, don't __dealloc__ here
        pass


    #_DataSet_2D(DataType tIn, int wIn, int pIn)

    #def virtual void Write3D(self,CpptrajFile, int, int, int) = 0 :

    #def virtual double GetElement(self,int, int, int) = 0 :

    #def virtual double operator[](self,size_t) = 0 :

    #def virtual size_t NX(self) = 0 :

    #def virtual size_t NY(self) = 0 :

    #def virtual size_t NZ(self) = 0 :

    #def void Add(self,size_t, void *):

    #def int Allocate_N_O_D(self,size_t, size_t, size_t, Vec3, Vec3):

    #def int Allocate_N_C_D(self,size_t, size_t, size_t, Vec3, Vec3):

    #def int Allocate_X_C_D(self, Vec3, Vec3, Vec3):

    #def int Allocate_N_O_Box(self,size_t, size_t, size_t, Vec3, Box):

    #def void GridInfo(self):

    #def bint CalcBins(self,double x, double y, double z, int i, int j, int k):

    #def void BinIndices(self,double x, double y, double z, int i, int j, int k):

    #def Vec3 BinCorner(self,int i, int j, int k):

    #def Vec3 BinCenter(self,int i, int j, int k):

    #def  Vec3 GridOrigin(self):

    #def Matrix_3x3 Ucell(self):

    #def double VoxelVolume(self):

