# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarray

cdef class DatasetMatrixDouble (DataSet_2D):
    def __cinit__(self):
        self.thisptr = new _DatasetMatrixDouble()
        self.baseptr_1 = <_DataSet_2D*> self.thisptr
        self.baseptr0 = <_DataSet*> self.thisptr

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def __getitem__(self, idx):
        return self.data[idx]

    def __iter__(self):
        for value in self.data:
            yield value

    def alloc(self):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = _DatasetMatrixDouble.Alloc()
        return dset

    @property
    def n_snapshots(self):
        return self.thisptr.Nsnapshots()

    def element(self, size_t x, size_t y):
        return self.thisptr.Element(x, y)

    def add_element(self, double d):
        return self.thisptr.AddElement(d)

    def set_element(self,size_t x, size_t y, double d):
        self.thisptr.SetElement(x, y, d)

    def vect(self):
        return self.thisptr.Vect()

    def allocate_vector(self,size_t vsize):
        self.thisptr.AllocateVector(vsize)

    def store_mass(self, Darray mIn):
        self.thisptr.StoreMass(mIn)

    @property
    def mass(self):
        return self.thisptr.Mass()

    def get_full_matrix(self):
        """return python array with length = n_rows*n_cols"""
        cdef int nr = self.n_rows
        cdef int nc = self.n_cols 
        cdef int i, j
        cdef pyarray arr0 = pyarray('d', [])

        for i in range(nr):
            for j in range(nc):
                arr0.append(self.baseptr_1.GetElement(i, j))
        return arr0

    @property
    def data(self):
        """return 1D python array of matrix' data"""
        return self.to_ndarray()

    def to_ndarray(self, copy=True):
        """use copy=True to be the same as DataSet_1D"""
        import numpy as np
        cdef int n_rows = self.n_rows
        cdef int n_cols = self.n_cols
        cdef double[:, :] dview = np.empty((n_rows, n_cols), dtype='f8')
        cdef int i, j

        for i in range(n_rows):
            for j in range(n_cols):
                dview[i, j] = self.baseptr_1.GetElement(i, j)
        return np.asarray(dview)

    def to_cpptraj_sparse_matrix(self):
        """return 1D numpy array, dtype='f8'
        """
        import numpy as np
        cdef int size = self.size
        cdef double[:] dview = np.empty(size, dtype='f8')

        for i in range(size):
            dview[i] = self.thisptr.index_opr(i)
        return np.asarray(dview)

    def to_half_matrix(self):
        import numpy as np
        hm = np.zeros((self.n_rows, self.n_cols)) 
        mt = self.to_cpptraj_sparse_matrix()

        hm[np.triu_indices(self.n_rows, 1)] = mt[mt !=0]
        return hm
