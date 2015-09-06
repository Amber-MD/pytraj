# distutils: language = c++


cdef class DataSet_3D (DataSet):
    def __cinit__(self):
        self.baseptr_1 = <_DataSet_3D*> self.baseptr0

    def __dealloc__(self):
        # since this is ABC, don't __dealloc__ here
        pass


    #_DataSet_2D(DataType tIn, int wIn, int pIn)

    #def virtual void Write3D(self,CpptrajFile, int, int, int) = 0 :

    #def virtual double GetElement(self,int, int, int) = 0 :

    #def virtual double operator[](self,size_t) = 0 :

    #def virtual size_t NX(self) = 0 :

    #def virtual size_t NY(self) = 0 :

    #def virtual size_t NZ(self) = 0 :

    #def void Add(self,size_t, void *):

    #def int Allocate_N_O_D(self,size_t, size_t, size_t, Vec3, Vec3):

    #def int Allocate_N_C_D(self,size_t, size_t, size_t, Vec3, Vec3):

    #def int Allocate_X_C_D(self, Vec3, Vec3, Vec3):

    #def int Allocate_N_O_Box(self,size_t, size_t, size_t, Vec3, Box):

    #def void GridInfo(self):

    #def bint CalcBins(self,double x, double y, double z, int i, int j, int k):

    #def void BinIndices(self,double x, double y, double z, int i, int j, int k):

    #def Vec3 BinCorner(self,int i, int j, int k):

    #def Vec3 BinCenter(self,int i, int j, int k):

    #def  Vec3 GridOrigin(self):

    #def Matrix_3x3 Ucell(self):

    #def double VoxelVolume(self):

# distutils: language = c++
from __future__ import absolute_import
from cython.view cimport array as cyarray
from ..utils import _import_numpy 

cdef class DatasetGridFloat(DataSet_3D):
    def __cinit__(self):
        self.baseptr0 = <_DataSet*> new _DatasetGridFloat()
        self.baseptr_1 = <_DataSet_3D*> self.baseptr0
        self.thisptr = <_DatasetGridFloat*> self.baseptr0
        self.py_free_mem = True

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def __str__(self):
        _, np = _import_numpy()
        basic_str = super(DataSet_3D, self).__str__() + "\n"
        if np:
            my_str = basic_str + "values: " + self.values.__str__()
        else:
            my_str = basic_str
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
        cdef float* ptr = &self.thisptr.index_opr(0)
        return <float[:nx, :ny, :nz]> ptr

    def to_ndarray(self, copy=True):
        # copy=True: is a dummy argument to be consistent with DataSet_1D
        has_np, np = _import_numpy()
        if not has_np:
            raise ImportError('require numpy')
        else:
            return np.array(self.data[:])

    def tolist(self):
        return [[list(x) for x in y] for y in self.data]
