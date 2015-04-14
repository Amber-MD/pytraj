# distutils: language = c++


cdef class Array1D:
    def __cinit__(self):
        self.thisptr = new _Array1D()

    def __dealloc__(self):
        del self.thisptr

    #def Array1D(self):

    #def Array1D(self, Array1D):

    #def Array1D(self, DataSetList):

    #def size_t DetermineMax(self):

    #def int push_back(self,DataSet_1D *):

    #def DataSet_1D * operator[](self,int idx):

    #def DataSet_1D * operator[](self,int idx):

    #def bint empty(self):

    #def _iterator begin(self):

    #def _iterator end(self):

    #def size_t size(self):

    #def void clear(self):

    #def void SortArray1D(self):

    #def int AddDataSets(self, DataSetList):

    #def int AddTorsionSets(self, DataSetList):

    #def int AddSetsFromArgs(self, ArgList, DataSetList):

    #def int CheckXDimension(self):

