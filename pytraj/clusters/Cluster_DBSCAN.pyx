# distutils: language = c++


cdef class Cluster_DBSCAN (ClusterList):
    def __cinit__(self):
        # self.thisptr: from ClusterList
        self.thisptr = new _Cluster_DBSCAN()
        self.ptr = <_Cluster_DBSCAN*>self.thisptr

    def __dealloc__(self):
        #del ptr
        del self.thisptr

    def help(self):
        self.ptr.Help()

    def setup_cluster(self,ArgList arglist):
        return self.ptr.SetupCluster(arglist.thisptr[0])

    def clustering_info(self):
        self.ptr.ClusteringInfo()

    def cluster(self):
        return self.ptr.Cluster()

    def add_sieved_frames(self):
        self.ptr.AddSievedFrames()

    def cluster_results(self,CpptrajFile cpptrajfile):
        self.ptr.ClusterResults(cpptrajfile.thisptr[0])

