# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from SymmetricRmsdCalc cimport *
from datasets.DataSet_Coords cimport *
from datasets.DataSet_1D cimport *
from clusters.ClusterMatrix cimport *

# TODO: double check declaration

ctypedef vector[int] Cframes
ctypedef vector[_DataSet*] DsArray

cdef extern from "ClusterDist.h": 
    #ctypedef vector[int] Cframes
    #ctypedef vector[_DataSet*] DsArray
    #ctypedef  Cframes "ClusterDist::Cframes" Cframes:
    #    pass
    cdef cppclass _Centroid "Centroid":
        pass
    cdef cppclass _Centroid_Coord "Centroid_Coord":
        _Centroid_Coord() 
        _Centroid_Coord(const _Frame& frame)
        _Centroid_Coord(int natom)
        _Centroid * Copy() 
    cdef cppclass _ClusterDist_Euclid "ClusterDist_Euclid":
        _ClusterDist_Euclid() 
        _ClusterDist_Euclid(const DsArray&)
        void PairwiseDist(_ClusterMatrix&, SievedFrames&)
        double FrameDist(int, int)
        double CentroidDist(_Centroid *, _Centroid *)
        double FrameCentroidDist(int, _Centroid *)
        void CalculateCentroid(_Centroid *, const Cframes&)
        _Centroid * NewCentroid(const Cframes&)
        _ClusterDist * Copy() 
    cdef cppclass _Centroid_Num "Centroid_Num":
        _Centroid_Num() 
        _Centroid_Num(double val)
        _Centroid * Copy() 
    cdef cppclass _ClusterDist_SRMSD "ClusterDist_SRMSD":
        _ClusterDist_SRMSD() 
        _ClusterDist_SRMSD(_DataSet *, const _AtomMask&, bint, bint, int)
        void PairwiseDist(_ClusterMatrix&, SievedFrames&)
        double FrameDist(int, int)
        double CentroidDist(_Centroid *, _Centroid *)
        double FrameCentroidDist(int, _Centroid *)
        void CalculateCentroid(_Centroid *, const Cframes&)
        _Centroid * New_Centroid(const Cframes&)
        _ClusterDist * Copy() 
    cdef cppclass _ClusterDist_Num "ClusterDist_Num":
        _ClusterDist_Num() 
        _ClusterDist_Num(_DataSet *)
        void PairwiseDist(_ClusterMatrix&, SievedFrames&)
        double FrameDist(int, int)
        double CentroidDist(_Centroid *, _Centroid *)
        double FrameCentroidDist(int, _Centroid *)
        void CalculateCentroid(_Centroid *, const Cframes&)
        _Centroid * NewCentroid(const Cframes&)
        _ClusterDist * Copy() 
    cdef cppclass _ClusterDist "ClusterDist":
        pass
    cdef cppclass _ClusterDist_DME "ClusterDist_DME":
        _ClusterDist_DME()
        _ClusterDist_DME(_DataSet *, const _AtomMask&)
        void PairwiseDist(_ClusterMatrix&, SievedFrames&)
        double FrameDist(int, int)
        double CentroidDist(_Centroid *, _Centroid *)
        double FrameCentroidDist(int, _Centroid *)
        void CalculateCentroid(_Centroid *, const Cframes&)
        _Centroid * NewCentroid(const Cframes&)
        _ClusterDist * Copy() 
    cdef cppclass _Centroid_Multi "Centroid_Multi":
        _Centroid_Multi() 
        _Centroid_Multi(const vector[double]& val)
        _Centroid * Copy() 
    cdef cppclass _ClusterDist_RMS "ClusterDist_RMS":
        _ClusterDist_RMS() 
        _ClusterDist_RMS(_DataSet *, const _AtomMask&, bint, bint)
        void PairwiseDist(_ClusterMatrix&, SievedFrames&)
        double FrameDist(int, int)
        double CentroidDist(_Centroid *, _Centroid *)
        double FrameCentroidDist(int, _Centroid *)
        void CalculateCentroid(_Centroid *, const Cframes&)
        _Centroid * NewCentroid(const Cframes&)
        _ClusterDist * Copy() 


cdef class Centroid_Coord:
    cdef _Centroid_Coord* thisptr

cdef class ClusterDist_Euclid:
    cdef _ClusterDist_Euclid* thisptr

cdef class Centroid_Num:
    cdef _Centroid_Num* thisptr

cdef class ClusterDist_SRMSD:
    cdef _ClusterDist_SRMSD* thisptr

cdef class Centroid:
    cdef _Centroid* thisptr

cdef class ClusterDist_Num:
    cdef _ClusterDist_Num* thisptr

cdef class ClusterDist:
    cdef _ClusterDist* thisptr

cdef class ClusterDist_DME:
    cdef _ClusterDist_DME* thisptr

cdef class Centroid_Multi:
    cdef _Centroid_Multi* thisptr

cdef class ClusterDist_RMS:
    cdef _ClusterDist_RMS* thisptr
