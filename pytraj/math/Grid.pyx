# distutils: language = c++
from __future__ import absolute_import
from cython.view cimport array as cyarray
from ..utils import _import_numpy 


cdef class Grid:
    def __cinit__(self, *args):
        cdef size_t x, y, z

        self.thisptr = new _Grid[float]()

        if len(args) == 3:
            x, y, z = args
            self.resize(x, y, z)

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, idx):
        cdef size_t x, y, z
        x, y, z = idx
        return self.thisptr.element(x, y, z)

    def __setitem__(self, idx, value):
        cdef size_t x, y, z
        x, y, z = idx
        self.thisptr.setGrid(x, y, z, <float> value)

    @property
    def size(self):
        return self.thisptr.size()

    def resize(self, size_t x, size_t y, size_t z):
        self.thisptr.resize(x, y, z)

    @property
    def nx(self):
        return self.thisptr.NX()

    @property
    def ny(self):
        return self.thisptr.NY()

    @property
    def nz(self):
        return self.thisptr.NZ()

    def _element(self, int x, int y, int z):
        return self.thisptr.element(x, y, z)

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
                    myview[i, j, k] = self.thisptr.element(i, j, k)
        return myview

    def to_ndarray(self):
        has_np, np = _import_numpy()
        if not has_np:
            raise ImportError("need numpy")
        else:
            return np.asarray(self.data[:], dtype=np.float32)

    def tolist(self):
        return [[list(x) for x in y] for y in self.data]
