# distutils: language = c++

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from libcpp.vector cimport vector
from cpython.array cimport array as pyarray
from ..utils.check_and_assert import _import_numpy, is_int, is_range
"""
In [1]: from Matrix_3x3 import Matrix_3x3 as M3x3
In [2]: import numpy as np
In [3]: x = np.arange(9).astype(float)
In [4]: m = M3x3(x)
In [5]: m.Print("3x3 matrix m: ")
    3x3 matrix m: 
       0.0000   1.0000   2.0000
       3.0000   4.0000   5.0000
       6.0000   7.0000   8.0000
In [6]: y = np.array([100.,]).astype(float)
In [7]: n = M3x3(y)
In [8]: n.Print("3x3 matrix n: ")
    3x3 matrix n: 
     100.0000 100.0000 100.0000
     100.0000 100.0000 100.0000
     100.0000 100.0000 100.0000
"""

cdef class Matrix_3x3:
    def __cinit__(self, Xin=None):
        """TODO: 
             add doc
             Add: mat1 = Matrix_3x3(mat2)
                  (Cython complains "TypeError: 'src.Matrix_3x3.Matrix_3x3' 
                   does not have the buffer interface")
        """
        cdef double[:] X 
        cdef double _xin = 0.
        import numbers
        import itertools

        if Xin is None:
            self.thisptr = new _Matrix_3x3(0.0)
        elif isinstance(Xin, numbers.Number):
            _xin = Xin
            self.thisptr = new _Matrix_3x3(_xin)
        else:
            if len(Xin) == 9:
                # 1D
                X = pyarray('d', Xin)
            elif len(Xin) == 3:
                # 2D: 3x3
                _flat_list = []
                if hasattr(Xin, 'tolist'):
                    Xin = Xin.tolist()
                for x in Xin:
                    for x0 in x:
                        _flat_list.append(x0)
                X = pyarray('d', _flat_list)
            if X.shape[0] == 9:
                #Takes array of 9, row-major
                self.thisptr = new _Matrix_3x3(&X[0])
            elif X.shape[0] == 1:
                #Set all elements to the same number
                self.thisptr = new _Matrix_3x3(X[0])
            elif X.shape[0] == 3:
                #Set Set diagonal
                x, y, z = X
                self.thisptr = new _Matrix_3x3(x, y, z)
            else: 
                raise ValueError("Must be array with length of None, 1, 3 or 9")

    def __iter__(self):
        cdef double[:] myview
        for myview in self.buffer2d[:]:
            yield myview

    def __getitem__(self, idx):
        return self.buffer2d[idx]

    def __setitem__(self, idx, value):
        self.buffer2d[idx] = value
    
    def __str__(self):
        txt = "Matrix 3x3: \n"
        for arr in self.buffer2d[:]:
            txt += " ".join(str(x) for x in arr)
            txt += "\n"
        return txt

    def __repr__(self):
        return self.__str__()

    def copy(self):
        cdef Matrix_3x3 other = Matrix_3x3()
        other.thisptr = new _Matrix_3x3(self.thisptr[0])
        return other

    def __dealloc__(self):
        """Free memory"""
        if self.thisptr:
            del self.thisptr

    def __imul__(Matrix_3x3 self, Matrix_3x3 other):
        """mat *= other"""

        self.thisptr[0].star_equal(other.thisptr[0])
        return self

    @property
    def row1(self):
        """
        Parameters: None
        ---------------

        Return
        ------
        Instance of Vec3
        """
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Row1()
        return vec

    @property
    def row2(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Row2()
        return vec

    @property
    def row3(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Row3()
        return vec

    @property
    def col1(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Col1()
        return vec

    @property
    def col2(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Col2()
        return vec

    @property
    def col3(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Col3()
        return vec

    def zeros(self):
        self.thisptr.Zero()

    def diagonalize(self, Vec3 vect): 
        self.thisptr.Diagonalize(vect.thisptr[0])

    def diagonalize_sort(self, Vec3 vectds):
        self.thisptr.Diagonalize_Sort(vectds.thisptr[0])

    def diagonalize_sort_chirality(self, Vec3 vectds, int idx):
        self.thisptr.Diagonalize_Sort_Chirality(vectds.thisptr[0], idx)

    def transpose(self):
        self.thisptr.Transpose()

    def rotation_around_zaxis(self, idx, idy):
        self.thisptr.RotationAroundZ(idx, idy)

    def rotation_around_yaxis(self, idx, idz):
        self.thisptr.RotationAroundY(idx, idz)

    def __mul__(Matrix_3x3 self, other):
        cdef Matrix_3x3 mat_other, mat_result
        cdef Vec3 vec_other, vec_result

        if isinstance(other, Vec3):
            vec_other = other
            vec_result = Vec3()
            vec_result.thisptr[0] = self.thisptr[0] * vec_other.thisptr[0]
            return vec_result
        elif isinstance(other, Matrix_3x3):
            mat_other = other
            mat_result = Matrix_3x3()
            mat_result.thisptr[0] = self.thisptr[0] * mat_other.thisptr[0]
            return mat_result
        else:
            raise ValueError("Must be either Matrix_3x3 or Vec3") 

    def calc_rotation_matrix(self, *args):
        """
        """
        cdef Vec3 vec
        cdef double theta
        cdef double x, y, z

        if len(args)== 2:
            vec, theta = args
            if not isinstance(vec, Vec3):
                raise ValueError("Must be a vector")
            self.thisptr.CalcRotationMatrix(vec.thisptr[0], theta)
        elif len(args) == 3:
            x, y, z = args
            self.thisptr.CalcRotationMatrix(x, y, z)
        else:
            raise ValueError('must be "Vec3, theta" or "x, y, z"')

    def rotation_angle(self):
       return self.thisptr.RotationAngle()

    def axis_of_rotation(self, theta):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.AxisOfRotation(theta)
        return vec 

    def transpose_mult(self, Matrix_3x3 other):
        cdef Matrix_3x3 result
        result =  Matrix_3x3()
        result.thisptr[0] = self.thisptr.TransposeMult(other.thisptr[0])
        return result 

    def tolist(self):
        return [list(x) for x in self.buffer2d[:]]

    def to_ndarray(self, copy=True):
        import numpy as np
        if copy:
            return np.array(self.buffer2d)
        else:
            return np.asarray(self.buffer2d)

    def as_ndmatrix(self):
        return self.to_ndmatrix()

    def to_ndmatrix(self):
        """convert to numpy matrix as a memory view. No data copy is made"""
        try:
            _, np = _import_numpy()
            mat = np.asmatrix(self.buffer2d[:])
            return mat
        except:
            raise ImportError("Must have numpy installed")

    property buffer2d:
        def __get__(self):
            cdef double[:, :] arr0 = <double[:3, :3]> self.thisptr.Dptr()
            return arr0

    property buffer1d:
        def __get__(self):
            cdef double[:] arr0 = <double[:9]> self.thisptr.Dptr()
            return arr0
