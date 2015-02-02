# distutils: language = c++


cdef class Centroid_Coord:
    def __cinit__(self, *args):
        cdef int natom
        cdef Frame frame 
        if not args:
            self.thisptr = new _Centroid_Coord()
        else:
            if len(args) == 1:
                if isinstance(args[0], Frame):
                    frame = args[0]
                    self.thisptr = new _Centroid_Coord(frame.thisptr[0])
                else:
                    natom = args[0]
                    self.thisptr = new _Centroid_Coord(natom)
            else:
                raise ValueError("")

    def __dealloc__(self):
        del self.thisptr

    def Copy(self):
        cdef Centroid cc = Centroid()
        cc.thisptr = self.thisptr.Copy()
        return cc

cdef class ClusterDist_Euclid:
    def __cinit__(self):
        self.thisptr = new _ClusterDist_Euclid()

    def __dealloc__(self):
        del self.thisptr

    #def ClusterDist_Euclid(self):

    #def ClusterDist_Euclid(self, DsArray):

    #def void PairwiseDist(self,ClusterMatrix, ClusterSieve:: SievedFrames):

    def FrameDist(self,int idx, int idy):
        return self.thisptr.FrameDist(idx, idy)

    #def double CentroidDist(self,Centroid *, Centroid *):

    #def double FrameCentroidDist(self,int, Centroid *):

    #def void CalculateCentroid(self,Centroid *, Cframes):

    #def Centroid * NewCentroid(self, Cframes):

    def Copy(self):
        # create a pointer to ClusterDist but don't create instance (abstract class)
        cdef ClusterDist cc  = ClusterDist()
        cc.thisptr = self.thisptr.Copy()
        return cc

cdef class Centroid_Num:
    def __cinit__(self, *args):
        pass
        #cdef double val
        #if not args:
        #    self.thisptr = new _Centroid_Num()
        #else:
        #    val = args[0]
        #    self.thisptr = new _Centroid_Num(val)

    def __dealloc__(self):
        del self.thisptr

    #def Copy(self):
    #    cdef Centroid cc = Centroid()
    #    cc.thisptr = self.thisptr.Copy()
    #    return cc

cdef class ClusterDist_SRMSD:
    def __cinit__(self):
        self.thisptr = new _ClusterDist_SRMSD()

    def __dealloc__(self):
        del self.thisptr

    #def ClusterDist_SRMSD(self):

    #def ClusterDist_SRMSD(self,DataSet *, AtomMask, bint, bint, int):

    #def void PairwiseDist(self,ClusterMatrix, ClusterSieve:: SievedFrames):

    #def double FrameDist(self,int, int):

    #def double CentroidDist(self,Centroid *, Centroid *):

    #def double FrameCentroidDist(self,int, Centroid *):

    #def void CalculateCentroid(self,Centroid *, Cframes):

    #def Centroid * NewCentroid(self, Cframes):

    #def ClusterDist * Copy(self):

cdef class Centroid:
    def __cinit__(self):
        pass
        #self.thisptr = new _Centroid()

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

cdef class ClusterDist_Num:
    def __cinit__(self):
        self.thisptr = new _ClusterDist_Num()

    def __dealloc__(self):
        del self.thisptr

    #def ClusterDist_Num(self):

    #def ClusterDist_Num(self,DataSet *):

    #def void PairwiseDist(self,ClusterMatrix, ClusterSieve:: SievedFrames):

    #def double FrameDist(self,int, int):

    #def double CentroidDist(self,Centroid *, Centroid *):

    #def double FrameCentroidDist(self,int, Centroid *):

    #def void CalculateCentroid(self,Centroid *, Cframes):

    #def Centroid * NewCentroid(self, Cframes):

    #def ClusterDist * Copy(self):

cdef class ClusterDist:
    def __cinit__(self):
        pass
        #self.thisptr = new _ClusterDist()

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

cdef class ClusterDist_DME:
    def __cinit__(self):
        self.thisptr = new _ClusterDist_DME()

    def __dealloc__(self):
        del self.thisptr

    #def ClusterDist_DME(self):

    #def ClusterDist_DME(self,DataSet *, AtomMask):

    #def void PairwiseDist(self,ClusterMatrix, ClusterSieve:: SievedFrames):

    #def double FrameDist(self,int, int):

    #def double CentroidDist(self,Centroid *, Centroid *):

    #def double FrameCentroidDist(self,int, Centroid *):

    #def void CalculateCentroid(self,Centroid *, Cframes):

    #def Centroid * NewCentroid(self, Cframes):

    #def ClusterDist * Copy(self):

cdef class Centroid_Multi:
    def __cinit__(self):
        self.thisptr = new _Centroid_Multi()

    def __dealloc__(self):
        del self.thisptr

    #def Centroid_Multi(self):

    #def Centroid_Multi(self, vector[double] val):

    #def Centroid * Copy(self):

cdef class ClusterDist_RMS:
    def __cinit__(self):
        self.thisptr = new _ClusterDist_RMS()

    def __dealloc__(self):
        del self.thisptr

    #def ClusterDist_RMS(self):

    #def ClusterDist_RMS(self,DataSet *, AtomMask, bint, bint):

    #def void PairwiseDist(self,ClusterMatrix, ClusterSieve:: SievedFrames):

    #def double FrameDist(self,int, int):

    #def double CentroidDist(self,Centroid *, Centroid *):

    #def double FrameCentroidDist(self,int, Centroid *):

    #def void CalculateCentroid(self,Centroid *, Cframes):

    #def Centroid * NewCentroid(self, Cframes):

    #def ClusterDist * Copy(self):

