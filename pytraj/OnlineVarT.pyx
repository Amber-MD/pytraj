# distutils: language = c++


cdef class StatsMap:
    def __cinit__(self):
        self.thisptr = new _StatsMap()

    def __dealloc__(self):
        del self.thisptr

    def accumulate(self,map a):
        self.thisptr.accumulate(a)

    def mean(self,int i):
        return self.thisptr.mean(i)

    def variance(self,int i):
        return self.thisptr.variance(i)

    def mean_generator(self):
        cdef iterator it = self.thisptr.mean_begin()
        while it != self.thisptr.mean_end():


    #def iterator mean_begin(self):

    #def iterator mean_end(self):

    #def iterator variance_begin(self):

    #def iterator variance_end(self):

    #def Value nData(self):

cdef class Stats:
    def __cinit__(self):
        self.thisptr = new _Stats()

    def __dealloc__(self):
        del self.thisptr

    #def Stats(self):

    #def void accumulate(self, Float x):

    #def Float mean(self):

    #def Float variance(self):

    #def Float nData(self):

