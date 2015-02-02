# distutils: language = c++
from ..ArgList cimport *


cdef extern from "ActionFrameCounter.h": 
    cdef cppclass _ActionFrameCounter "ActionFrameCounter":
        _Action_FrameCounter() 
        int InitFrameCounter(_ArgList&)
        bint CheckFrameCounter(int frameNum) const 
        void _FrameCounterInfo() const 
        void _FrameCounterBrief() const 


cdef class ActionFrameCounter:
    cdef _ActionFrameCounter* thisptr

