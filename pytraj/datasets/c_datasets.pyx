# distutils: language = c++

from __future__ import division
from ..c_dict import DataTypeDict, ScalarTypeDict, get_key
from ..datafiles.datafiles import DataFileList, DataFile

import operator
from cython.operator cimport preincrement as incr, dereference as deref
from cython.view cimport array as cyarray

# python level
import numpy as np
from ..utils import is_int
from ..shared_methods import _xyz, my_str_method
from ..cyutils import get_positive_idx
from ..c_traj.c_trajectory import TrajectoryCpptraj
from ..topology cimport Topology
from ..externals.six import string_types


cdef class Dataset:
    __cpptraj_dataset__ = None
    """
    Original doc from cpptraj
    -------------------------
    Class: Dataset
        Base class that all Dataset types will inherit.
        Datasets are given certain attributes to make Dataset selection easier;
        these are name, index, and aspect. Name is typically associated with the
        action that creates the Dataset, e.g. RMSD or distance. Index is used
        when and action outputs subsets of data, e.g. with RMSD it is possible to
        output per-residue RMSD, where the Dataset index corresponds to the residue
        number. Aspect is used to further subdivide output data type; e.g. with
        nucleic acid analysis each base pair (divided by index) has shear,
        stagger etc calculated.

    pytraj doc
    ----------
        Add some other methods: `tolist`, `to_ndarray`, `plot` and attribute `data`

    """

    def __cinit__(self):
        pass
        # don't create instance for this abstract class
        # self.baseptr0 = new _Dataset()

    def __dealloc__(self):
        pass
        # let sub-class do this job
        # if self.baseptr0 != NULL:
        #    del self.baseptr0

    property name:
        def __get__(self):
            name = self.baseptr0.Meta().Name()
            return name.decode()

        def __set__(self, name):
            cdef _MetaData meta = self.baseptr0.Meta()
            meta.SetName(name.encode())
            self.baseptr0.SetMeta(meta)

    property aspect:
        def __get__(self):
            aspect = self.baseptr0.Meta().Aspect()
            return aspect.decode()

        def __set__(self, aspect):
            cdef string s = aspect.encode()
            cdef _MetaData meta = self.baseptr0.Meta()

            meta.SetAspect(s)
            self.baseptr0.SetMeta(meta)

    property _legend:
        def __get__(self):
            _legend = self.baseptr0.Meta().Legend()
            return _legend.decode()

        def __set__(self, _legend):
            cdef string s = _legend.encode()
            self.baseptr0.SetLegend(s)

    property key:
        # retire self._legend?
        def __get__(self):
            return self._legend

        def __set__(self, _legend):
            self._legend = _legend

    property dtype:
        def __get__(self):
            return get_key(self.baseptr0.Type(), DataTypeDict).lower()

    def __str__(self):
        cname = self._class_name
        size = self.size
        _legend = self._legend
        msg0 = """<pytraj.datasets.{0}: size={1}, key={2}> """.format(
            cname, size, _legend)
        return msg0

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        raise NotImplementedError("Must over-write Dataset data attr")

    def __getitem__(self, idx):
        raise NotImplementedError("Must over-write Dataset data attr")

    def __setitem__(self, idx, value):
        raise NotImplementedError("Must over-write Dataset data attr")

    def __array__(self):
        """
        Aim: directly use numpy to perform analysis without casting to ndararay again

        Examples
        -------
            d = Dataset_integer()
            d.resize(200)
            d.data[:] = np.arange(200)
            np.mean(d)
        """
        try:
            return np.asarray(self.data)
        except:
            raise NotImplementedError("don't know how to cast to ndarray")

    @property
    def _class_name(self):
        return self.__class__.__name__

    @property
    def size(self):
        return self.baseptr0.Size()

    @property
    def data(self):
        """mostly memoryview
        """
        raise NotImplementedError("Must over-write Dataset data attr")

    property scalar_type:
        def __set__(self, stype):
            '''
            '''
            cdef _MetaData meta = self.baseptr0.Meta()
            cdef scalarType st

            try:
                st = ScalarTypeDict[stype.upper()]
            except KeyError:
                raise KeyError(ScalarTypeDict.keys())

            meta.SetScalarType(ScalarTypeDict[stype.upper()])
            self.baseptr0.SetMeta(meta)

        def __get__(self):
            return get_key(self.baseptr0.Meta().ScalarType(), ScalarTypeDict)

    def tolist(self):
        return list(self.data)

    property values:
        '''return a copy or non-copy, depending on data
        '''

        def __set__(self, values):
            self.data = values

        def __get__(self):
            return np.asarray(self.data)

    def to_ndarray(self, copy=False):
        """return ndarray view of self.data"""
        if not copy:
            return np.asarray(self.data)
        else:
            # make copy
            return np.array(self.data)

    def to_dict(self, use_numpy=False):
        if np and use_numpy:
            return {self._legend: self.values}
        if not np and use_numpy:
            raise ImportError("require numpy. Set `use_numpy=False`")
        return {self._legend: self.tolist()}


cdef class Dataset1D (Dataset):
    def __cinit__(self, *args):
        cdef Dataset dset
        # make sure two pointers pointing to the same address
        self.baseptr_1 = <_Dataset1D*> self.baseptr0

    def __dealloc__(self):
        pass

    def __str__(self):
        basic_str = super(Dataset1D, self).__str__() + "\n"
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
            self.baseptr_1 = <_Dataset1D*> self.baseptr0
        elif idx == 1:
            self.baseptr0 = <_Dataset*> self.baseptr_1
        else:
            raise ValueError("idx must be 0 or 1")

    def _allocate_1D(self, size_t size):
        cdef vector[size_t] v
        v.push_back(size)
        return self.baseptr_1.Allocate(v)

    def from_array_like(self, array_like):
        """
        Notes: require numpy
        """
        old_size = self.size
        self.resize(self.size + len(array_like))
        self.values[old_size:] = array_like

    def _xcrd(self):
        '''x-data.
        '''
        cdef unsigned int idx

        return np.array([self.baseptr_1.Xcrd(idx) for idx in range(len(self))])


cdef class DatasetDouble (Dataset1D):
    def __cinit__(self, *args):
        # TODO : Use only one pointer?
        self.baseptr0 = <_Dataset*> new _DatasetDouble()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_Dataset1D*> self.baseptr0
        self.thisptr = <_DatasetDouble*> self.baseptr0

        # let Python/Cython free memory
        self._own_memory = True

        if args:
            if isinstance(args[0], list):
                self.data = args[0]

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def __getitem__(self, idx):
        # return self.thisptr.index_opr(idx)
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

    def append(self, double d):
        self.thisptr.AddElement(d)

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

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

cdef class DatasetFloat (Dataset1D):
    def __cinit__(self):
        self.baseptr0 = <_Dataset*> new _DatasetFloat()
        self.baseptr_1 = <_Dataset1D*> self.baseptr0
        self.thisptr = <_DatasetFloat*> self.baseptr0

        # let Python/Cython free memory
        self._own_memory = True

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

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
            cdef float x
            self.resize(0)
            for x in data:
                self.thisptr.AddElement(x)

    def append(self, ds):
        cdef int new_size = self.size + ds.size
        cdef int j
        self.resize(new_size)

        j = 0
        for i in range(self.size, new_size):
            self[i] = ds[j]
            j += 1


cdef class DatasetInteger (Dataset1D):
    def __cinit__(self):
        # TODO : Use only one pointer?
        self.baseptr0 = <_Dataset*> new _DatasetInteger()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_Dataset1D*> self.baseptr0
        self.thisptr = <_DatasetInteger*> self.baseptr0

        # let Python/Cython free memory
        self._own_memory = True

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def __getitem__(self, idx):
        # return self.thisptr.index_opr(idx)
        cdef int i

        if is_int(idx):
            return self.thisptr.index_opr(idx)
        elif isinstance(idx, slice):
            if idx == slice(None):
                return np.array([self.thisptr.index_opr(i)
                                 for i in range(self.size)])
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

    def append(self, val):
        self.thisptr.AddElement(<int> val)

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
            cdef size_t size = len(data)
            cdef unsigned int i

            self.thisptr.Resize(size)
            # let numpy handle, just need to resize self
            values = np.asarray(self.data)
            values[:] = data


cdef class DatasetString (Dataset1D):
    def __cinit__(self):
        self.baseptr0 = <_Dataset*> new _DatasetString()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_Dataset1D*> self.baseptr0
        self.thisptr = <_DatasetString*> self.baseptr0

        # let Python/Cython free memory
        self._own_memory = True

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

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

    property data:
        def __get__(self):
            return [s.decode() for s in self]

        def __set__(self, data):
            self.thisptr.Resize(len(data))
            for i, x in enumerate(data):
                # x must be string
                self[i] = x.encode()

    def tolist(self):
        return self.data


cdef class DatasetVector(Dataset):
    def __cinit__(self):
        self._own_memory = True
        self.thisptr = new _DatasetVector()
        self.baseptr0 = <_Dataset*> self.thisptr

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    @property
    def shape(self):
        return (self.size, 3)

    def __getitem__(self, idx):
        """return memoryview for Vec3. No data is copied.
        """
        cdef Vec3 vec = Vec3()
        if idx == -1:
            idx = self.size - 1
        vec._own_memory = False
        vec.thisptr = &(self.thisptr.index_opr(idx))
        return vec

    def __iter__(self):
        for i in range(self.size):
            yield self[i]

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

    def append(self, Vec3 vec):
        self.thisptr.AddVxyz(vec.thisptr[0])

    property data:
        def __get__(self):
            cdef int i
            cdef int size = self.size
            cdef _Vec3 _vec3
            cdef double[:, ::1] dview = np.empty((size, 3), dtype='f8')

            # copy data to arr by dview
            for i in range(size):
                _vec3 = self.thisptr.index_opr(i)
                dview[i, 0] = _vec3.Dptr()[0]
                dview[i, 1] = _vec3.Dptr()[1]
                dview[i, 2] = _vec3.Dptr()[2]
            return np.asarray(dview, dtype='f8')

        def __set__(self, values):
            '''values must be 2D array that support memory view (such as numpy array)
            '''
            cdef int i
            cdef double[:] xyz
            cdef _Vec3 _vec
            cdef double[:, ::1] arr = values

            if arr.shape[1] != 3:
                raise ValueError("must have shape = (n_frames, 3))")

            self.resize(0)
            for i in range(arr.shape[0]):
                xyz = arr[i]
                _vec.Assign( &xyz[0])
                self.thisptr.AddVxyz(_vec)

    def tolist(self):
        # overwrite
        # x is memview array
        return [x.tolist() for x in self.data]

    def to_ndarray(self, copy=True):
        return np.asarray(self.data)

    def to_dataframe(self):
        import pandas as pd
        return pd.DataFrame(self.to_ndarray(), columns=list('xyz'))


cdef class Dataset2D (Dataset):
    def __cinit__(self):
        # since Dataset2D inherits from Dataset, make sure two pointers pointing
        # to the same address
        self.baseptr_1 = <_Dataset2D*> self.baseptr0

    def __dealloc__(self):
        pass

    property kind:
        def __get__(Dataset2D self):
            '''
            '''
            cdef int i = <int> self.baseptr_1.Kind()

            kind_dict = {0: 'full', 1: 'half', 2: 'tri'}
            return kind_dict[i]

    @property
    def n_rows(self):
        return self.baseptr_1.Nrows()

    @property
    def n_cols(self):
        return self.baseptr_1.Ncols()

    def get_element(self, int x, int y):
        return self.baseptr_1.GetElement(x, y)

    def _allocate_2D(self, size_t x, size_t y):
        cdef vector[size_t] v
        v.push_back(x)
        v.push_back(y)
        self.baseptr_1.Allocate(v)

    def _allocate_half(self, size_t x):
        self.baseptr_1.AllocateHalf(x)

    def _allocate_triangle(self, size_t x):
        self.baseptr_1.AllocateTriangle(x)

    def get_full_matrix(self):
        raise NotImplementedError("must over-write in subclass")

    def to_dataframe(self):
        raise NotImplementedError("must overwrite in subclass")

cdef class DatasetMatrixDouble (Dataset2D):
    def __cinit__(self):
        self.thisptr = new _DatasetMatrixDouble()
        self.baseptr_1 = <_Dataset2D*> self.thisptr
        self.baseptr0 = <_Dataset*> self.thisptr

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def __getitem__(self, idx):
        return self.data[idx]

    def __iter__(self):
        for value in self.data:
            yield value

    @property
    def n_snapshots(self):
        return self.thisptr.Nsnapshots()

    def element(self, size_t x, size_t y):
        return self.thisptr.Element(x, y)

    def append(self, double d):
        return self.thisptr.AddElement(d)

    def set_element(self, size_t x, size_t y, double d):
        self.thisptr.SetElement(x, y, d)

    def vect(self):
        return self.thisptr.Vect()

    def _allocate_vector(self, size_t vsize):
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
        cdef double[:, ::1] arr = np.empty((nr, nc), dtype='f8')

        for i in range(nr):
            for j in range(nc):
                arr[i, j] = self.baseptr_1.GetElement(i, j)
        return np.asarray(arr)

    property data:
        def __get__(self):
            """return 1D python array of matrix' data"""
            return self.to_ndarray()

    def _set_data_half_matrix(self, double[:] values, size_t vsize, size_t n_cols):
        '''only support half matrix
        TODO: correct?
        '''
        cdef unsigned int i

        (<_Dataset2D*> self.thisptr).AllocateHalf(n_cols)

        for i in range(vsize):
            self.thisptr.AddElement(values[i])

    def to_ndarray(self, copy=True):
        """use copy=True to be the same as Dataset1D"""
        cdef int n_rows = self.n_rows
        cdef int n_cols = self.n_cols
        cdef double[:, :] dview = np.empty((n_rows, n_cols), dtype='f8')
        cdef int i, j

        for i in range(n_rows):
            for j in range(n_cols):
                dview[i, j] = self.baseptr_1.GetElement(i, j)
        return np.asarray(dview)

    def _to_cpptraj_sparse_matrix(self):
        """return 1D numpy array, dtype='f8'
        """
        cdef int size = self.size
        cdef double[:] dview = np.empty(size, dtype='f8')

        for i in range(size):
            dview[i] = self.thisptr.index_opr(i)
        return np.array(dview)

    def to_half_matrix(self):
        hm = np.zeros((self.n_rows, self.n_cols))
        mt = self._to_cpptraj_sparse_matrix()

        hm[np.triu_indices(self.n_rows, 1)] = mt[mt !=0]
        return hm

cdef class DatasetMatrixFloat (Dataset2D):
    def __cinit__(self):
        self.thisptr = new _DatasetMatrixFloat()
        self.baseptr_1 = <_Dataset2D*> self.thisptr
        self.baseptr0 = <_Dataset*> self.thisptr

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def __getitem__(self, idx):
        return self.data[idx]

    def get_full_matrix(self):
        cdef int nr = self.n_rows
        cdef int nc = self.n_cols
        cdef int i, j
        cdef float[:, ::1] arr = np.empty((nr, nc), dtype='f4')

        for i in range(nr):
            for j in range(nc):
                arr[i, j] = self.baseptr_1.GetElement(i, j)
        return np.asarray(arr)

    @property
    def data(self):
        return self.get_full_matrix()

    def to_ndarray(self, copy=True):
        # use copy=True to be consistent with Dataset1D
        arr = np.array(self.get_full_matrix()).reshape(
            self.n_rows, self.n_cols)
        return arr

    def to_ndarray(self, copy=True):
        """use copy=True to be the same as Dataset1D"""
        cdef int n_rows = self.n_rows
        cdef int n_cols = self.n_cols
        cdef float[:, :] dview = np.empty((n_rows, n_cols), dtype='f4')
        cdef int i, j

        for i in range(n_rows):
            for j in range(n_cols):
                dview[i, j] = self.baseptr_1.GetElement(i, j)
        return np.asarray(dview)

    def _to_cpptraj_sparse_matrix(self):
        """return 1D numpy array, dtype='f8'
        """
        cdef int size = self.size
        cdef float[:] dview = np.empty(size, dtype='f4')

        for i in range(size):
            dview[i] = self.thisptr.index_opr(i)
        return np.asarray(dview)

    def to_half_matrix(self):
        hm = np.zeros((self.n_rows, self.n_cols))
        mt = self._to_cpptraj_sparse_matrix()

        hm[np.triu_indices(self.n_rows, 1)] = mt[mt !=0]
        return hm


cdef class Dataset3D (Dataset):
    def __cinit__(self):
        self.baseptr_1 = <_Dataset3D*> self.baseptr0

    def __dealloc__(self):
        # since this is ABC, don't __dealloc__ here
        pass

cdef class DatasetGridFloat(Dataset3D):
    def __cinit__(self):
        self.baseptr0 = <_Dataset*> new _DatasetGridFloat()
        self.baseptr_1 = <_Dataset3D*> self.baseptr0
        self.thisptr = <_DatasetGridFloat*> self.baseptr0
        self._own_memory = True

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def __str__(self):
        basic_str = super(Dataset3D, self).__str__() + "\n"
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

    property data:
        def __get__(self):
            """return a copy of 3D array of Grid"""
            cdef size_t nx, ny, nz
            nx, ny, nz = self.nx, self.ny, self.nz
            cdef float* ptr = &self.thisptr.index_opr(0)
            return <float[:nx, :ny, :nz]> ptr

        def __set__(self, float[:, :, :] values):
            cdef unsigned int nx, ny, nz
            cdef unsigned int i, j, k

            # use for old cython
            nx, ny, nz = [_ for _ in values.shape[:3]]
            self.resize(nx, ny, nz)

            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        self.thisptr.SetElement(i, j, k, values[i, j, k])

    def to_ndarray(self, copy=True):
        # copy=True: is a dummy argument to be consistent with Dataset1D
        return np.array(self.data[:])

    def tolist(self):
        return [[list(x) for x in y] for y in self.data]

cdef class DatasetModes(Dataset):
    def __cinit__(self):
        self._own_memory = True
        self.thisptr = new _DatasetModes()
        self.baseptr0 = <_Dataset*> self.thisptr

    def __dealloc__(self):
        if self.thisptr and self._own_memory:
            del self.thisptr

    property data:
        def __get__(self):
            return self.eigenvalues, self.eigenvectors

    property values:
        def __get__(self):
            return self.data

    property n_modes:
        def __get__(self):
            return self.thisptr.Nmodes()

    def eigval_to_freq(self, x):
        return self.thisptr.EigvalToFreq(<double> x)

    @property
    def vector_size(self):
        return self.thisptr.VectorSize()

    def _is_reduced(self):
        return self.thisptr.IsReduced()

    def _set_modes(self, bint is_reduced, int n_modes, int vsize, double[:] eigenvalues, double[:] eigenvectors):
        '''
        Notes
        -----
        eigenvectors is 1D array, make sure to reshape if yours is 2D
        '''
        self.thisptr.SetModes(is_reduced, n_modes, vsize, &eigenvalues[0], &eigenvectors[0])

    property eigenvalues:
        def __get__(self):
            cdef int i
            return np.array([self.thisptr.Eigenvalue(i) for i in
                             range(self.thisptr.Nmodes())])

    property eigenvectors:
        def __get__(self):
            cdef const double * ptr = self.thisptr.Eigenvectors()
            cdef int n_modes = self.thisptr.Nmodes()
            cdef int vsize = self.vector_size

            return np.array([ptr[i] for i in
                             range(n_modes*vsize)]).reshape(n_modes, vsize)

    def _allocate_avgcoords(self, int n):
        self.thisptr.AllocateAvgCoords(n)

    def _set_avg_frame(self, double[:] arr):
        cdef unsigned int i 
        cdef double* ptr = self.thisptr.AvgFramePtr() 

        for i in range(arr.shape[0]):
            ptr[i] = arr[i]

    def _get_avg_crd(self):
        return self.thisptr.AvgCrd()

cdef class DatasetMatrix3x3 (Dataset):
    def __cinit__(self):
        # TODO : Use only one pointer?
        self.baseptr0 = <_Dataset*> new _DatasetMatrix3x3()
        # make sure 3 pointers pointing to the same address?
        self.thisptr = <_DatasetMatrix3x3*> self.baseptr0

        # let Python/Cython free memory
        self._own_memory = True

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def __getitem__(self, int idx):
        '''return a copy
        '''
        if self.size <= 0:
            raise ValueError("size should be > 0")

        cdef Matrix_3x3 mat = Matrix_3x3()
        mat.thisptr[0] = self.thisptr[0][idx]
        return mat

    def __setitem__(self, int idx, double value):
        raise NotImplementedError('does not support setitem')

    def __iter__(self):
        """return copy"""
        if self.size <= 0:
            raise ValueError("size should be > 0")
        cdef vector[_Matrix_3x3].iterator it = self.thisptr.begin()
        cdef Matrix_3x3 mat

        while it != self.thisptr.end():
            mat = Matrix_3x3()
            mat.thisptr[0] = deref(it)
            incr(it)
            yield mat

    def append(self, Matrix_3x3 mat):
        self.thisptr.AddMat3x3(mat.thisptr[0])

    def tolist(self):
        return self.to_ndarray().tolist()

    property data:
        def __get__(self):
            return np.array([np.array(x) for x in self])

    def to_ndarray(self, copy=True):
        """return a copy
        """
        return np.asarray(self.data)

    def _append_from_array(self, arr):
        '''arr can be 2D array, shape=(n_frames, 9) or 3D array, shape=(n_frames, 3, 3)
        '''
        values = np.asarray(arr)
        if values.ndim == 2:
            self._append_from_2array(values)
        elif values.ndim == 3:
            self._append_from_2array(
                values.reshape(
                    values.shape[0],
                    values.shape[1] *
                    values.shape[2]))

    def _append_from_2array(self, double[:, ::1] arr):
        cdef unsigned int i

        for i in range(arr.shape[0]):
            self.thisptr.AddMat3x3(_Matrix_3x3( &arr[i, 0]))

cdef class DatasetMesh (Dataset1D):
    def __cinit__(self):
        self.baseptr0 = <_Dataset*> new _DatasetMesh()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_Dataset1D*> self.baseptr0
        self.thisptr = <_DatasetMesh*> self.baseptr0

        # let Python/Cython free memory
        self._own_memory = True

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def tolist(self):
        """return 2D list with format [index, value]
        """
        cdef unsigned int i
        return [[self.thisptr.X(i), self.thisptr.Y(i)] for i in range(self.size)]

    property data:
        def __get__(self):
            cdef double[:, ::1] dview = np.empty((self.size, 2), dtype='f8')
            cdef unsigned int i
            cdef unsigned int size = self.size

            # fill data for arr by using its dview
            for i in range(size):
                dview[i, 0], dview[i, 1] = self.thisptr.X(i), self.thisptr.Y(i)
            return np.asarray(dview)

    def _append_from_array(self, double[:, ::1] values):
        cdef unsigned int i

        for i in range(values.shape[0]):
            self.thisptr.AddXY(values[i, 0], values[i, 1])

    def to_ndarray(self, copy=True):
        """use copy=True to make consistent with Dataset1D
        """
        return np.asarray(self.data)

cdef class DatasetCoords(Dataset):
    def __cinit__(self):
        # abstract class, dont' create new object here
        # pass
        # make sure that two pointers pointing to the same address
        self.baseptr0 = <_Dataset*> self.baseptr_1
        self._top = Topology()

    def __dealloc__(self):
        # abstract class
        pass

    @property
    def n_frames(self):
        return self.size

    @property
    def n_atoms(self):
        """used for frame_iter"""
        return self.top.n_atoms

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        """iterately getting Frame instance
        TODO : get memoryview or copy?
        """
        cdef unsigned int i
        cdef unsigned int size = self.size
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.baseptr_1.AllocateFrame()

        for i in range(size):
            self.baseptr_1.GetFrame(i, frame.thisptr[0])
            yield frame

    def __getitem__(self, idx):
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.baseptr_1.AllocateFrame()

        if self.size == 0:
            raise ValueError("Your Trajectory is empty, how can I index it?")
        self.baseptr_1.GetFrame(idx, frame.thisptr[0])
        self.tmpfarray = frame
        return self.tmpfarray

    def __setitem__(self, int idx, Frame other):
        idx_1 = get_positive_idx(idx, self.size)
        # raise index out of range
        if idx != 0 and idx_1 == 0:
            # need to check if array has only 1 element.
            # arr[0] is  arr[-1]
            if idx != -1:
                raise ValueError("index is out of range")
        self.baseptr_1.SetCRD(idx, other.thisptr[0])

    def allocate_frame(self):
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.baseptr_1.AllocateFrame()
        return frame

    property top:
        def __get__(self):
            self._top.thisptr[0] = self.baseptr_1.Top()
            return self._top

        def __set__(self, Topology other):
            # self.baseptr_1.SetTopology(other.thisptr[0])
            self.baseptr_1.CoordsSetup(other.thisptr[0], self.baseptr_1.CoordsInfo())

    def add_frame(self, Frame frame):
        if self.top.n_atoms != frame.n_atoms:
            raise ValueError("Frame and Topology must have the same number of atoms")
        self.baseptr_1.AddFrame(frame.thisptr[0])

    def append(self, frame):
        """alis of addframe"""
        self.add_frame(frame)

    def get_frame(self, int idx, Frame frameout):
        self.baseptr_1.GetFrame(idx, frameout.thisptr[0])

    @property
    def xyz(self):
        """return a copy of xyz coordinates (ndarray, shape=(n_frames, n_atoms, 3)
        We can not return a memoryview since Trajectory is a C++ vector of Frame object
        """
        cdef Frame frame
        cdef int i, n_frames, n_atoms
        cdef double[:, :, ::1] xyz

        n_frames = self.n_frames
        n_atoms = self.top.n_atoms
        xyz = np.empty((n_frames, n_atoms, 3), dtype='f8')

        frame = Frame(n_atoms, xyz[0], _as_ptr=True)

        for i in range(n_frames):
            # use `frame` as a pointer pointing to `xyz` memory
            # dump coords to xyz array
            frame.thisptr.SetXptr(n_atoms, &xyz[i, 0, 0])
            # copy coordinates of `self[i]` to j-th frame in `traj`
            self.baseptr_1.GetFrame(i, frame.thisptr[0])
        return np.asarray(xyz)

    def tolist(self):
        """return flatten list for traj-like object"""
        cdef Frame frame
        return [frame.tolist() for frame in self]

    def to_dataframe(self):
        raise NotImplementedError()

cdef class DatasetCoordsCRD (DatasetCoords):
    def __cinit__(self):
        self.thisptr = new _DatasetCoordsCRD()
        self.baseptr0 = <_Dataset*> self.thisptr
        self.baseptr_1 = <_DatasetCoords*> self.thisptr

        # let python frees memory
        self._own_memory = True

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def load(self, filename):
        trajin = TrajectoryCpptraj()
        trajin._load(filename, self.top)

        for frame in trajin:
            self.append(frame.copy())

    property data:
        def __get__(self):
            return self

    property values:
        def __get__(self):
            return self.data

cdef class DatasetCoordsRef (DatasetCoords):
    def __cinit__(self):
        self.thisptr = new _DatasetCoordsRef()
        self.baseptr0 = <_Dataset*> self.thisptr
        self.baseptr_1 = <_DatasetCoords*> self.thisptr

        # let python frees memory
        self._own_memory = True

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def get_frame(self):
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.thisptr.RefFrame()
        # use self._base to avoid eager memory-free
        self._base = frame
        return self._base

    @property
    def values(self):
        """"""
        return self.data.xyz

    property data:
        def __get__(self):
            """"""
            return self.get_frame()

        def __set__(self, Frame frame):
            self.thisptr.SetCRD(0, frame.thisptr[0])

    property xyz:
        def __get__(self):
            # use np.array to make a copy to avoid memory free
            return np.array(self.data.xyz)


cdef class DatasetTopology(Dataset):
    def __cinit__(self):
        self.thisptr = new _DatasetTopology()
        self.baseptr0 = <_Dataset*> self.thisptr

        # let python frees memory
        self._own_memory = True

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    property _top:
        def __get__(self):
            cdef Topology top = Topology()
            top.thisptr[0] = self.thisptr.Top()
            return top

        def __set__(self, Topology top):
            self.thisptr.SetTop(top.thisptr[0])

    @property
    def values(self):
        """"""
        return self.data

    property data:
        def __get__(self):
            """"""
            return self._top

        def __set__(self, Topology top):
            self._top = top
