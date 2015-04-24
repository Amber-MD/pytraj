# distutils: language = c++

cdef class Dimension:
    def __cinit__(self, *args):
        """Add more args to function"""
        # create temp variables
        cdef double dmin, dmax
        cdef int dstep, dbins
        cdef string label
        cdef Dimension rhs

        if not args:
            self.thisptr = new _Dimension()
        elif len(args) == 1:
            # copy
            rhs = args[0]
            self.thisptr = new _Dimension(rhs.thisptr[0])
        elif len(args) == 3:
            dmin, dstep, dbins = args
            self.thisptr = new _Dimension(dmin, dstep, dbins)
        elif len(args) == 4:
            dmin, dstep, dbins, label = args
            self.thisptr = new _Dimension(dmin, dstep, dbins, label)
        else:
            raise ValueError("Not implemented")

    def __dealloc__(self):
        del self.thisptr

    property Label:
        def __get__(self):
            return self.thisptr.Label()
        def __set__(self, string label):
            self.thisptr.SetLabel(label)

    property Min:
        def __get__(self):
            return self.thisptr.Min()
        def __set__(self, double x):
            self.thisptr.SetMin(x)

    property Max:
        def __get__(self):
            return self.thisptr.Max()
        def __set__(self, x):
            self.thisptr.SetMax(x)
    
    property Step:
        def __get__(self):
            return self.thisptr.Step()
        def __set__(self, double x):
            self.thisptr.SetStep(x)

    property Bins:
        def __get__(self):
            return self.thisptr.Bins()
        def __set__(self, int bins):
            self.thisptr.SetBins(bins)

    def MinIsSet(self):
        return self.thisptr.MinIsSet()

    def MaxIsSet(self):
        return self.thisptr.MaxIsSet()

    def Coord(self, size_t x):
        return self.thisptr.Coord(x)

    def CalcBinsOrStep(self):
        return self.thisptr.CalcBinsOrStep()

    def PrintDim(self):
        self.thisptr.PrintDim()
    
