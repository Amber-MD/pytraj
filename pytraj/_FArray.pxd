# distutils: language = c++
from libcpp.vector cimport vector
from pytraj.Frame cimport _Frame

ctypedef vector[_Frame].iterator _FArray_iter
cdef extern from "FrameArray.h": 
    cdef cppclass _FArray "FrameArray":
        #_FrameArray() 
        #void resize(int nIn)
        _Frame& operator[](int idx)
        #void append "AddFrame"(const _Frame& fIn)
        #int SetupFrames(const vector[_Atom]& _Atoms, const _CoordinateInfo& cInfo)
        _FArray_iter begin() 
        _FArray_iter end() 
