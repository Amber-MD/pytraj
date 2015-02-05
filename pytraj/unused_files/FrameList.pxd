# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.TopologyList cimport _TopologyList, TopologyList
from pytraj.ReferenceFrame cimport _ReferenceFrame, ReferenceFrame
from pytraj.Frame cimport _Frame, Frame


cdef extern from "FrameList.h": 
    cdef cppclass _FrameList "FrameList":
        _FrameList() 
        #~_FrameList() 
        void Clear() 
        void SetDebug(int)
        _Frame ActiveReference() const 
        int SetActiveRef(int)
        int AddRefFrame(_ArgList&, const _TopologyList&)
        _ReferenceFrame GetFrameFromArgs(_ArgList&) const 
        _ReferenceFrame GetFrameByName(const string&) const 
        void List() const 
        int NumFrames() const 
        const char* RefArgs

cdef class FrameList:
    cdef _FrameList* thisptr
    #cdef const char* info
    cdef bint py_free_mem
