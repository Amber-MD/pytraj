# distutils: language = c++
from __future__ import absolute_import
from cython.view cimport array as cyarray
from ..utils import _import_numpy 
from ..exceptions import PytrajError

cdef class DataSet_GridFlt(DataSet_3D):
    def __cinit__(self):
        self.baseptr0 = <_DataSet*> new _DataSet_GridFlt()
        self.baseptr_1 = <_DataSet_3D*> self.baseptr0
        self.thisptr = <_DataSet_GridFlt*> self.baseptr0
        self.py_free_mem = True

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def __str__(self):
        _, np = _import_numpy()
        basic_str = super(DataSet_3D, self).__str__() + "\n"
        if np:
            try:
                my_str = basic_str + "values: " + self.values.__str__()
            except:
                my_str = basic_str
        else:
            my_str = basic_str + "(install numpy for pretty print)"
        return my_str

    def __getitem__(self, idx):
        cdef size_t x, y, z
        x, y, z = idx
        return self.thisptr.GetElement(x, y, z)

    def __setitem__(self, idx, value):
        cdef size_t x, y, z
        x, y, z = idx
        self.thisptr.SetElement(x, y, z, <float> value)

    def resize(self, size_t x, size_t y, size_t z):
        self.thisptr.Allocate3D(x, y, z)

    @property
    def nx(self):
        return self.thisptr.NX()

    @property
    def ny(self):
        return self.thisptr.NY()

    @property
    def nz(self):
        return self.thisptr.NZ()

    @property
    def shape(self):
        return (self.nx, self.ny, self.nz)

    @property
    def data(self):
        """return a copy of 3D array of Grid"""
        cdef size_t nx, ny, nz
        nx, ny, nz = self.nx, self.ny, self.nz
        cdef int i, j, k
        cdef cyarray carr = cyarray(shape=(nx, ny, nz),
                                   itemsize=sizeof(float), format="f")
        cdef float[:, :, :] myview = carr

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    myview[i, j, k] = self.thisptr.GetElement(i, j, k)
        return myview

    def to_ndarray(self):
        has_np, np = _import_numpy()
        if not has_np:
            raise PytrajError("need numpy")
        else:
            return np.array(self.data[:])

    def tolist(self):
        return [[list(x) for x in y] for y in self.data]
