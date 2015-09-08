# distutils: language = c++
from cpython.array cimport array as pyarray
from cython.view cimport array as cyarray
from ..utils import _import_numpy

cdef class DatasetMesh (DataSet_1D):
    def __cinit__(self):
        self.baseptr0 = <_DataSet*> new _DatasetMesh()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.thisptr = <_DatasetMesh*> self.baseptr0

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

    def tolist(self):
        """return 2D list with format [index, value]
        """
        # xcrd is for cpptraj's output which use index starting of 1
        # we need to subtract "1"
        return [[int(self._xcrd(i)-1), self._d_val(i)] for i in range(self.size)]

    def to_ndarray(self, copy=True):
        """use copy=True to make consistent with DataSet_1D
        """
        _, np = _import_numpy()
        return np.array(self.tolist())
