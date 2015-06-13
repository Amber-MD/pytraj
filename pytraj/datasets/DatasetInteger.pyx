# distutils: language = c++
from cpython.array cimport array as pyarray
from cython.view cimport array as cyarray

# python level
from ..utils import is_int

cdef class DatasetInteger (DataSet_1D):
    def __cinit__(self):
        # TODO : Use only one pointer? 
        self.baseptr0 = <_DataSet*> new _DatasetInteger()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.thisptr = <_DatasetInteger*> self.baseptr0

        # let Python/Cython free memory
        self.py_free_mem = True

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
        cdef pyarray arr0 = pyarray('i', [])
        cdef int i

        if is_int(idx):
            return self.thisptr.index_opr(idx)
        elif isinstance(idx, slice):
            if idx == slice(None):
                for i in range(self.size):
                    arr0.append(self.thisptr.index_opr(i))
                return arr0
            else:
                raise NotImplementedError("only support slice(None)")
        else:
            raise NotImplementedError("only support single indexing or slice(None)")

    def __setitem__(self, int idx, int value):
        cdef int * ptr
        ptr = &(self.thisptr.index_opr(idx))
        ptr[0] = value
        
    def __iter__(self):
        cdef int i
        cdef int size = self.size
        for i in range(size):
            yield self.thisptr.index_opr(i)

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

    def count(self, value=None):
        """
        Parameters
        value : int, optional

        Examples
        --------
        ds.count()
        ds.count(1)
        """
        cdef int i, count

        if value is None:
            from collections import Counter
            return Counter(self.data)
        else:
            count = 0
            for i in self:
                if value == i:
                    count += 1
            return count

    def append(self, values):
        cdef int i, d
        cdef int[:] int_view
        cdef pyarray arr

        if hasattr(values, 'real') and hasattr(values, 'imag'):
            # a number
            self.thisptr.AddElement(<int> values)
        else:
            try:
                int_view = values
            except:
                if hasattr(values, 'data'):
                    try:
                        int_view = values.data
                    except:
                        arr = pyarray('i', values)
                        int_view = arr

            for i in range(int_view.shape[0]):
                self.thisptr.AddElement(int_view[i])

    def _add(self, int idx, int value):
        self.thisptr.Add(idx, &value)

    property data:
        def __get__(self):
            """return memoryview of data array
            """
            cdef cyarray myview
            cdef int size = self.size
            cdef int* ptr

            if size == 0:
                return None
            ptr = &self.thisptr.index_opr(0)
            myview = <int[:size]> ptr
            return myview

        def __set__(self, data):
            cdef vector[int] v
            cdef int x
            cdef size_t size = len(data)

            self.baseptr_1.Allocate1D(size)
            self.data[:] = data
