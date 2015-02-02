# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.Topology cimport _Topology, Topology
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.AtomMask cimport _AtomMask, AtomMask
from pytraj.FileName cimport _FileName, FileName
from pytraj.Frame cimport _Frame, Frame


cdef extern from "ReferenceFrame.h": 
    cdef cppclass _ReferenceFrame "ReferenceFrame":
        _ReferenceFrame() 
        _ReferenceFrame(int)
        #~_ReferenceFrame() 
        const _Frame& Coord() const 
        const _Topology& Parm() const 
        bint error() const 
        bint empty() const 
        const _FileName& FrameName() const 
        const string& Tag() const 
        int LoadRef(const string&, _Topology *, int)
        int LoadRef(const string&, _ArgList&, _Topology *, const string&, int)
        int StripRef(const _AtomMask&)
        void RefInfo() const 
        void ClearRef() 


cdef class ReferenceFrame:
    cdef _ReferenceFrame* thisptr

