# distutils: language = c++


cdef class DataSet_Modes (DataSet):
    def __cinit__(self):
        self.thisptr = new _DataSet_Modes()
        self.baseptr0 = <_DataSet*> self.thisptr

    def __dealloc__(self):
        del self.thisptr

    #def DataSet_Modes(self):

    #def DataSet * Alloc(self):

    #def size_t Size(self):

    #def int Sync(self):

    #def void Info(self):

    #def void Add(self,size_t, void *):

    #def AvgIt AvgBegin(self):

    #def  Darray AvgCrd(self):

    #def  Darray Mass(self):

    #def int NavgCrd(self):

    #def double * AvgFramePtr(self):

    #def  double * AvgFramePtr(self):

    #def void AllocateAvgCoords(self,int n):

    #def void SetAvgCoords(self, DataSet_2D):

    #def int SetModes(self,bint, int, int, double *, double *):

    #def int CalcEigen(self, DataSet_2D, int):

    #def void PrintModes(self):

    #def int EigvalToFreq(self,double):

    #def int MassWtEigvect(self,DataSet_MatrixDbl:: Darray):

    #def int ReduceVectors(self):

    #def int Thermo(self,CpptrajFile, int, double, double):

    #def void SetType(self,DataSet_2D::MatrixType typeIn):

    #def double Eigenvalue(self,int i):

    #def  double * Eigenvectors(self):

    #def  double * Eigenvector(self,int i):

    #def int Nmodes(self):

    #def int VectorSize(self):

    #def DataSet_2D::MatrixType Type(self):

    #def bint IsReduced(self):

