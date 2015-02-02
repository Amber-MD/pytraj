# distutils: language = c++


cdef class ActionFrameCounter:
    def __cinit__(self):
        self.thisptr = new _ActionFrameCounter()

    def __dealloc__(self):
        del self.thisptr

    #def int InitFrameCounter(self,ArgList):

    #def bint CheckFrameCounter(self,int frameNum):

    #def void FrameCounterInfo(self):

    #def void FrameCounterBrief(self):

