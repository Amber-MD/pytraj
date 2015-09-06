# distutils: language = c++
from __future__ import division
from pytraj.utils import _import_numpy
import operator


cdef class DataSet_1D (DataSet):
    def __cinit__(self, *args):
        cdef DataSet dset
        # make sure two pointers pointing to the same address
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0

    def __dealloc__(self):
        pass

    def __str__(self):
        _, np = _import_numpy()
        basic_str = super(DataSet_1D, self).__str__() + "\n"
        if np:
            my_str = basic_str + "values: \n" + self.values.__str__()
        else:
            my_str = basic_str
        return my_str

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.size

    @property
    def shape(self):
        return (self.size,)

    def _recast_pointers(self, idx=0):
        """
        Since we use >=2 pointers pointing to the same address,
        we need to recast after each pointer assignment
        """
        if idx == 0:
            self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        elif idx == 1:
            self.baseptr0 = <_DataSet*> self.baseptr_1
        else:
            raise ValueError("idx must be 0 or 1")

    def allocate_1D(self, size_t sizet):
        return self.baseptr_1.Allocate1D(sizet)

    def _d_val(self, size_t sizet):
        return self.baseptr_1.Dval(sizet)

    def _xcrd(self, size_t sizet):
        return self.baseptr_1.Xcrd(sizet)

    def _is_torsion_array(self):
        return self.baseptr_1.IsTorsionArray()

    def from_array_like(self, array_like):
        """
        Notes: require numpy
        """
        old_size = self.size
        self.resize(self.size + len(array_like))
        self.values[old_size:] = array_like

    def avg(self):
        return sum(self.values) / len(self)

    def mean(self, *args, **kwd):
        import numpy as np
        return np.mean(self.values, *args, **kwd)

    def mean_with_error(self, DataSet other):
        m0 = self.mean()
        m1 = other.mean() 
        return ((m0 + m1)/2., abs(m0 - m1)/2.)

    def min(self, *args, **kwd):
        return self.baseptr_1.Min()

    def max(self, *args, **kwd):
        return self.baseptr_1.Max()

    def cross_corr(self, DataSet_1D D2, DataSet_1D Ct, int lagmaxIn, 
                bint calccovar, bint usefft):
        return self.baseptr_1.CrossCorr(D2.baseptr_1[0], Ct.baseptr_1[0], 
                lagmaxIn, calccovar, usefft)

    def corr_coeff(self, DataSet_1D other):
        return self.baseptr_1.CorrCoeff(other.baseptr_1[0])

    # below are copied from `dask` package: New BSD
    # see pytraj/licenses/externals/dask.txt for license
    def __abs__(self):
        return elemwise(operator.abs, self)
    def __add__(self, other):
        return elemwise(operator.add, self, other)
    def __radd__(self, other):
        return elemwise(operator.add, other, self)
    def __and__(self, other):
        return elemwise(operator.and_, self, other)
    def __rand__(self, other):
        return elemwise(operator.and_, other, self)
    def __div__(self, other):
        return elemwise(operator.div, self, other)
    def __rdiv__(self, other):
        return elemwise(operator.div, other, self)
    def __invert__(self):
        return elemwise(operator.invert, self)
    def __lshift__(self, other):
        return elemwise(operator.lshift, self, other)
    def __rlshift__(self, other):
        return elemwise(operator.lshift, other, self)
    def __mod__(self, other):
        return elemwise(operator.mod, self, other)
    def __rmod__(self, other):
        return elemwise(operator.mod, other, self)
    def __mul__(self, other):
        return elemwise(operator.mul, self, other)
    def __rmul__(self, other):
        return elemwise(operator.mul, other, self)
    def __neg__(self):
        return elemwise(operator.neg, self)
    def __or__(self, other):
        return elemwise(operator.or_, self, other)
    def __pos__(self):
        return self
    def __ror__(self, other):
        return elemwise(operator.or_, other, self)
    def __rpow__(self, other):
        return elemwise(operator.pow, other, self)
    def __rshift__(self, other):
        return elemwise(operator.rshift, self, other)
    def __rrshift__(self, other):
        return elemwise(operator.rshift, other, self)
    def __sub__(self, other):
        return elemwise(operator.sub, self, other)
    def __rsub__(self, other):
        return elemwise(operator.sub, other, self)
    def __truediv__(self, other):
        return elemwise(operator.truediv, self, other)
    def __rtruediv__(self, other):
        return elemwise(operator.truediv, other, self)
    def __floordiv__(self, other):
        return elemwise(operator.floordiv, self, other)
    def __rfloordiv__(self, other):
        return elemwise(operator.floordiv, other, self)
    def __xor__(self, other):
        return elemwise(operator.xor, self, other)
    def __rxor__(self, other):
        return elemwise(operator.xor, other, self)

    # end of copy from dask

def elemwise(op, self, other=None):
    if other:
        if hasattr(other, 'values'):
            _other = other.values
        else:
            _other = other
        if hasattr(self, 'values'):
            _self = self.values
        else:
            _self = self
        return op(_self, _other)
    else:
        return op(self.values)
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
# distutils: language = c++
from cpython.array cimport array as pyarray
from cython.view cimport array as cyarray


cdef class DatasetFloat (DataSet_1D):
    def __cinit__(self):
        self.baseptr0 = <_DataSet*> new _DatasetFloat()
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.thisptr = <_DatasetFloat*> self.baseptr0

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
        return self.data[idx]

    def __setitem__(self, idx, value):
        self.data[idx] = value
        
    def __iter__(self):
        cdef int i
        for i in range(self.size):
            yield self.thisptr.index_opr(i)

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

    property data:
        def __get__(self):
            """return memoryview of data array
            """
            cdef cyarray myview
            cdef int size = self.size
            cdef float* ptr

            if size == 0:
                return None
            ptr = &self.thisptr.index_opr(0)
            myview = <float[:size]> ptr
            return myview

        def __set__(self, data):
            raise NotImplementedError()

    def append(self, ds):
        cdef int new_size = self.size + ds.size
        cdef int j
        self.resize(new_size)

        j = 0
        for i in range(self.size, new_size):
            self[i] = ds[j]
            j += 1
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
# distutils: language = c++
from cpython.array cimport array as pyarray

# python level
#from pytraj.optional_libs import HAS_NUMPY, ndarray

cdef class DatasetString (DataSet_1D):
    def __cinit__(self):
        self.baseptr0 = <_DataSet*> new _DatasetString()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.thisptr = <_DatasetString*> self.baseptr0

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

    def __getitem__(self, int idx):
        return self.thisptr.index_opr(idx)

    def __setitem__(self, int idx, value):
        cdef string* ptr
        ptr = &(self.thisptr.index_opr(idx))
        ptr[0] = value
        
    def __iter__(self):
        cdef int i
        for i in range(self.size):
            yield self.thisptr.index_opr(i)

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

    @property
    def data(self):
        return [s.decode() for s in self]

    def tolist(self):
        return self.data

    def to_pyarray(self):
        cdef pyarray arr0 = pyarray('u', self.tolist())
        return arr0
