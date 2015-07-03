# distutils: language = c++
from cpython.array cimport array as pyarray
from cython.view cimport array as cyarray

cdef class DatasetDouble (DataSet_1D):
    def __cinit__(self, *args):
        # TODO : Use only one pointer? 
        self.baseptr0 = <_DataSet*> new _DatasetDouble()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.thisptr = <_DatasetDouble*> self.baseptr0

        # let Python/Cython free memory
        self.py_free_mem = True

        if args:
            if isinstance(args[0], list):
                self.data = args[0]

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr


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

    def add_element(self, double d):
        self.thisptr.AddElement(d)

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

    def info(self):
        self.thisptr.Info()

    def xcrd(self, size_t idx):
        raise NotImplementedError()

    def append(self, dset, idx=None):
        cdef DatasetDouble dset_
        cdef double elm
        cdef size_t idx_

        if isinstance(dset, DatasetDouble):
            if idx is not None:
                raise ValueError("can not use id with DatasetDouble instance")
            dset_ = dset
            self.thisptr.Append(dset_.thisptr[0])
        else:
            # try to add a `double` elm
            elm = dset
            idx_ = <size_t> idx
            self.thisptr.Add(idx_, <void*> (&elm))

    property data:
        def __get__(self):
            """return memoryview of data array
            """
            cdef cyarray myview
            cdef int size = self.size
            cdef double* ptr

            if size == 0:
                return None
            ptr = &self.thisptr.index_opr(0)
            myview = <double[:size]> ptr
            return myview

        def __set__(self, data):
            cdef vector[double] v
            cdef double x

            for x in data:
                # really need to do this?
                v.push_back(<double> x)
            self.thisptr.assign_opr(v)
