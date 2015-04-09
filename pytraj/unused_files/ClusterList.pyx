# distutils: language = c++


cdef class ClusterList:
    def __cinit__(self):
        # Don't create new instance for this class
        # It has virtual methods
        pass

    def __dealloc__(self):
        # Let the sub-class do its own destructor?
        #del self.thisptr
        pass

    def metric_string(self, DistMetricType dm):
        return self.thisptr.MetricString(dm)

    def nclusters(self):
        return self.thisptr.Nclusters()

    def set_debug(self,int debugIn):
        self.thisptr.SetDebug(debugIn)

    def renumber(self,bint addSievedFrames):
        self.thisptr.Renumber(addSievedFrames)

    def summary(self, string summaryfile, int maxframesIn):
        self.thisptr.Summary(summaryfile, maxframesIn)

    def summary_part(self, string summaryfile, int maxframesIn, vector[int] splitFrames):
        # TODO: 
        self.thisptr.Summary_Part(summaryfile, maxframesIn, splitFrames)

    def print_clusters_to_file(self, string filename, int maxframesIn):
        self.thisptr.PrintClustersToFile(filename, maxframesIn)

    def print_clusters(self):
        self.thisptr.PrintClusters()

    # need to wrap DsArray
    #def calc_frame_distances(self, string filename, DsArray dataSetslist, DistModeType mode, DistMetricType metric, bint nofit, bint useMass, string maskexpr, int sieve, int sieveSeed):
    #    return self.thisptr.CalcFrameDistances(filename, dataSets, mode, metric, nofit, useMass, maskexpr, sieve, sieveSeed)

    #def  cluster_iterator begincluster(self):
    #def  cluster_iterator endcluster(self):

