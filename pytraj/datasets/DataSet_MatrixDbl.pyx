# distutils: language = c++


cdef class DataSet_MatrixDbl (DataSet_2D):
    def __cinit__(self):
        self.thisptr = new _DataSet_MatrixDbl()
        self.baseptr_1 = <_DataSet_2D*> self.thisptr
        self.baseptr0 = <_DataSet*> self.thisptr

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def __getitem__(self, int idx):
        if idx >= self.size:
            # double-check
            # what is self.size?
            raise IndexError("out of index")
        return self.thisptr.index_opr(idx)

    def alloc(self):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = _DataSet_MatrixDbl.Alloc()
        return dset

    @property
    def size(self):
        return self.thisptr.Size()

    def info(self):
        self.thisptr.Info()

    #def allocate2D(self,size_t x, size_t y):
    #    return self.thisptr.Allocate2D(x, y)

    #def allocate_half(self,size_t x):
    #    return self.thisptr.AllocateHalf(x)

    #def allocate_triangle(self, size_t x):
    #    return self.thisptr.AllocateTriangle(x)

    #def write2D(self, CpptrajFile cppfile, int dix, int idy):
    #    self.thisptr.Write2D(cppfile, idx, idy)

    #def get_element(self,size_t x, size_t y):
    #    return self.thisptr.GetElement(x, y)

    #def double * MatrixArray(self):
    #def MatrixKind Kind(self):
    #def MatrixType Type(self):

    def element(self, size_t x, size_t y):
        return self.thisptr.Element(x, y)

    def add_element(self, double d):
        return self.thisptr.AddElement(d)

    def set_element(self,size_t x, size_t y, double d):
        self.thisptr.SetElement(x, y, d)

    #def iterator begin(self):
    #def iterator end(self):

    def vect(self):
        return self.thisptr.Vect()

    #def v1(self):
    #    return self.thisptr.V1()

    def allocate_vector(self,size_t vsize):
        self.thisptr.AllocateVector(vsize)

    #def Darray::iterator v1begin(self):
    #def Darray::iterator v1end(self):
    #def void SetTypeAndKind(self,MatrixType tIn, MatrixKind kIn):

    def store_mass(self, Darray mIn):
        self.thisptr.StoreMass(mIn)

    def mass(self):
        return self.thisptr.Mass()
