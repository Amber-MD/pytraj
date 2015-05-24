# distutils: language = c++
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
            my_str = basic_str
        return my_str

    def __repr__(self):
        return self.__str__()

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

    def write_buffer(self, CpptrajFile cppfile, size_t sizet):
        self.baseptr_1.WriteBuffer(cppfile.thisptr[0], sizet)

    def _d_val(self, size_t sizet):
        return self.baseptr_1.Dval(sizet)

    def _xcrd(self, size_t sizet):
        return self.baseptr_1.Xcrd(sizet)

    def _is_torsion_array(self):
        return self.baseptr_1.IsTorsionArray()

    def avg(self, *args):
        if not args:
            return self.baseptr_1.Avg()
        else:
            sd = args[0]
            return self.baseptr_1.Avg(sd)

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
