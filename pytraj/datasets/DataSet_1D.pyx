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
