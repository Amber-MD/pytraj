# distutils: language = c++
from __future__ import division
from pytraj.utils import _import_numpy


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
            try:
                my_str = basic_str + "values: " + self.values.__str__()
            except:
                my_str = basic_str
        else:
            my_str = basic_str + "(install numpy for pretty print)"
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

    def mean(self):
        return self.avg()

    def mean_with_error(self, DataSet other):
        m0 = self.mean()
        m1 = other.mean() 
        return ((m0 + m1)/2., abs(m0 - m1)/2.)

    def min(self):
        return self.baseptr_1.Min()

    def max(self):
        return self.baseptr_1.Max()

    def cross_corr(self, DataSet_1D D2, DataSet_1D Ct, int lagmaxIn, 
                bint calccovar, bint usefft):
        return self.baseptr_1.CrossCorr(D2.baseptr_1[0], Ct.baseptr_1[0], 
                lagmaxIn, calccovar, usefft)

    def corr_coeff(self, DataSet_1D other):
        return self.baseptr_1.CorrCoeff(other.baseptr_1[0])
