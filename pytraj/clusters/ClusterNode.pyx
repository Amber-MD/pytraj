# distutils: language = c++


cdef class ClusterNode:
    def __cinit__(self):
        self.thisptr = new _ClusterNode()

    def __dealloc__(self):
        del self.thisptr

    #def ClusterNode(self):

    #def ClusterNode(self,ClusterDist *, ClusterDist:: Cframes, int):

    #def ClusterNode(self, ClusterNode):

    #def  bint operator[(self, ClusterNode):

    #def  void MergeFrames(self, ClusterNode):

    #def int FindBestRepFrame(self, ClusterMatrix):

    #def void CalcEccentricity(self, ClusterMatrix):

    #def void CalculateCentroid(self,ClusterDist * Cdist):

    #def double CalcAvgToCentroid(self,ClusterDist *):

    #def frame_iterator beginframe(self):

    #def frame_iterator endframe(self):

    #def int ClusterFrame(self,int idx):

    #def  double AvgDist(self):

    #def  double Eccentricity(self):

    #def  int Num(self):

    #def  int Nframes(self):

    #def  int BestRepFrame(self):

    #def  Centroid * Cent(self):

    #def void SetAvgDist(self,double avg):

    #def void AddFrameToCluster(self,int fnum):

    #def void SetNum(self,int numIn):

    #def void SortFrameList(self):

    #def bint HasFrame(self,int):

    #def void RemoveFrameFromCluster(self,int):

    #def void RemoveFrameUpdateCentroid(self,ClusterDist *, int):

    #def void AddFrameUpdateCentroid(self,ClusterDist *, int):

