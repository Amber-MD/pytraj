# distutils: language = c++

from cython.view cimport array as cyarray
from pytraj.utils import _import
from pytraj._xyz import XYZ


cdef class DataSet_Vector (DataSet_1D):
    def __cinit__(self):
        self.py_free_mem = True
        self.thisptr = new _DataSet_Vector()
        self.baseptr0 = <_DataSet*> self.thisptr
        self.baseptr_1= <_DataSet_1D*> self.thisptr

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def __str__(self):
        _, pd = _import("pandas")
        if pd:
            return self.to_dataframe().__str__()
        else:
            print ("don't have pandas: use simple __str__")
            return super(DataSet_Vector, self).__str__()

    def alloc(self):
        cdef DataSet d0 = DataSet()
        d0.baseptr0 = self.thisptr.Alloc()
        return d0

    def __getitem__(self, idx):
        """return memoryview for Vec3. No data is copied.
        """
        cdef Vec3 vec = Vec3()
        if idx == -1:
            idx = self.size - 1
        vec.py_free_mem = False
        vec.thisptr = &(self.thisptr.index_opr(idx))
        return vec

    def __iter__(self):
        for i in range (self.size):
            yield self[i]

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

    def append(self, Vec3 vec):
        self.thisptr.AddVxyz(vec.thisptr[0])

    def tolist(self):
        # overwrite
        # x is memview array
        return [x.tolist() for x in self.data]

    def to_ndarray(self):
        import numpy as np
        # overwrite
        # x is memview array
        return np.asarray([np.asarray(x[:]) for x in self.data])

    def to_dataframe(self):
        from pytraj.utils import _import
        _, pd = _import("pandas")
        if pd:
            return pd.DataFrame(self.to_ndarray(), columns=list('xyz'))

    @property
    def data(self):
        """return self.__iter__
        Not sure what else we should return
        """
        return self.__iter__()

    def is_ired(self):
        return self.thisptr.IsIred()

    def set_ired(self):
        self.thisptr.SetIred()

    @property
    def values(self):
        return XYZ(self.to_ndarray())
