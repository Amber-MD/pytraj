# distutils: language = c++
from cpython.array cimport array as pyarray
from pytraj.cpptraj_dict import MatrixDict, MatrixKindDict, get_key

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

    @property
    def mtype(self):
        return get_key(self.thisptr.matType(), MatrixDict)

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
