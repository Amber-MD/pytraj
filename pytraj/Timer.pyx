# distutils: language = c++


cdef class Timer:
    def __cinit__(self):
        self.thisptr = new _Timer()

    def __dealloc__(self):
        del self.thisptr

    def start(self):
        self.thisptr.Start()

    def stop(self):
        self.thisptr.Stop()

    def total(self):
        return self.thisptr.Total()

    #def WriteTiming(self,int, char *, double):

    #def WriteTiming(self,int i, char * h):

