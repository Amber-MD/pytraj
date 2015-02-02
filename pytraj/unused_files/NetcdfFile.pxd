# distutils: language = c++
from libcpp.string cimport string
from pycpptraj.ReplicaDimArray cimport *


cdef extern from "NetcdfFile.h": 
#cdef extern from "NetcdfFile.cpp": 
    # NetcdfFile.h
    ctypedef enum NCTYPE "NetcdfFile::NCTYPE":
        NC_UNKNOWN "NetcdfFile::NC_UNKNOWN"
        NC_AMBERTRAJ "NetcdfFile::NC_AMBERTRAJ"
        NC_AMBERRESTART "NetcdfFile::NC_AMBERRESTART"
    cdef cppclass _NetcdfFile "NetcdfFile":
        NCTYPE GetNetcdfConventions(const char *)
        _NetcdfFile() 
        void NetcdfDebug() 
        string GetAttrText(const char *)
        #NCTYPE GetNetcdfConventions() 
        int NC_openRead(const string&)
        int NC_openWrite(const string&)
        int NC_createReservoir(bint, double, int, int&, int&)
        int NC_create(const string&, NCTYPE, int, bint, bint, bint, bint, bint, bint, const _ReplicaDimArray&, const string&)
        void NC_close() 
        int SetupFrame() 
        int SetupCoordsVelo(bint)
        int SetupTime() 
        int SetupBox(double *, NCTYPE)
        int SetupTemperature() 
        int SetupMultiD(_ReplicaDimArray&)
        void FloatToDouble(double *, const float *)
        void DoubleToFloat(float *, const double *)
        inline int Ncid() 
        inline int Ncatom() 
        inline int Ncatom3() 
        inline int Ncframe() 
        bint HasVelocities() 
        bint HasCoords() 
        inline void SetNcatom(int natomIn)


cdef class NetcdfFile:
    cdef _NetcdfFile* thisptr

