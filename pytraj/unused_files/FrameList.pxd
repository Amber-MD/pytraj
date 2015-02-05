# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.TopologyList cimport _TopologyList, TopologyList
from pytraj.ReferenceFrame cimport _ReferenceFrame, ReferenceFrame
from pytraj.Frame cimport _Frame, Frame

cdef cppclass _FrameList:
    _FrameList()
    __dealloc__()

cdef class FrameList:
    cdef _FrameList* thisptr
