# distutils: language = c++
from cpython.array cimport array as pyarray

# python level
#from pytraj.optional_libs import HAS_NUMPY, ndarray

cdef class DataSet_double (DataSet_1D):
    def __cinit__(self, *args):
        # TODO : Use only one pointer? 
        self.baseptr0 = <_DataSet*> new _DataSet_double()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.thisptr = <_DataSet_double*> self.baseptr0

        # let Python/Cython free memory
        self.py_free_mem = True

        if args:
            if isinstance(args[0], list):
                self.data = args[0]

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    #def DataSet_double(self):

    #@classmethod
    #def alloc(cls):
    #    '''return a memoryview as DataSet instane'''
    #    cdef DataSet dset = DataSet()
    #    dset.baseptr0 = _DataSet_double.Alloc()
    #    return dset

    def alloc(self):
        '''return a memoryview as DataSet instane'''
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.Alloc()
        return dset

    def __getitem__(self, idx):
        #return self.thisptr.index_opr(idx)
        # use self.data so we can use fancy indexing
        return self.data[idx]

    def __setitem__(self, int idx, double value):
        cdef double* ptr
        ptr = &(self.thisptr.index_opr(idx))
        ptr[0] = value
        
    def __iter__(self):
        cdef int i
        for i in range(self.size):
            yield self.thisptr.index_opr(i)

    property data:
        def __get__(self):
            """return a copy of data
            # (return a Python list, 
            Cython convert vector[double]) to list[double]
            """
            return self.thisptr.Data()

        def __set__(self, data):
            cdef vector[double] v
            cdef double x

            for x in data:
                # really need to do this?
                v.push_back(<double> x)
            self.thisptr.assign_opr(v)

    def add_element(self, double d):
        self.thisptr.AddElement(d)

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

    # turn of this property to use base class' method
    #@property
    #def size(self):
    #    return self.thisptr._Size()

    def info(self):
        self.thisptr.Info()

    # already declare in base-class
    #def allocate1D(self, size_t size):
    #    self.thisptr.Allocate1D(size)

    #def void Add(self,size_t, void *):
    #def double Dval(self,size_t idx):

    def xcrd(self, size_t idx):
        raise NotImplementedError()

    # already declare in base-class
    #def write_buffer(self, CpptrajFile cpptrajfile, size_t idx):
    #    self.thisptr.WriteBuffer(cpptrajfile.thisptr[0], idx)

    def append(self, dset, idx=None):
        cdef DataSet_double dset_
        cdef double elm
        cdef size_t idx_

        if isinstance(dset, DataSet_double):
            if idx is not None:
                raise ValueError("can not use id with DataSet_double instance")
            dset_ = dset
            self.thisptr.Append(dset_.thisptr[0])
        else:
            # try to add a `double` elm
            elm = dset
            idx_ = <size_t> idx
            self.thisptr.Add(idx_, <void*> (&elm))

    #def void SetNOE(self,double b, double bh, double r):
    #def double NOE_bound(self):
    #def double NOE_boundH(self):
    #def double NOE_rexp(self):
    #def void ShiftTorsions(self,double, double):
