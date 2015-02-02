# distutils: language = c++

cdef class AnalysisList:
    def __cinit__(self):
        self.thisptr = new _AnalysisList()

    def __dealloc__(self):
        del self.thisptr

    #def AnalysisList(self):

    #def void Clear(self):

    #def void SetDebug(self,int):

    #def int Debug(self):

    #def int AddAnalysis(self,DispatchObject::DispatchAllocatorType, ArgList, TopologyList *, DataSetList *, DataFileList *):

    #def int DoAnalyses(self):

    #def void List(self):

    #def bint Empty(self):

