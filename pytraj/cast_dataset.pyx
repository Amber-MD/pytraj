# distutils: language = c++
from pytraj.datasets.DataSet_1D cimport DataSet_1D, _DataSet_1D
from pytraj.datasets.DataSet_2D cimport DataSet_2D, _DataSet_2D
from pytraj.datasets.DataSet_double cimport DataSet_double, _DataSet_double
from pytraj.datasets.DataSet_float cimport DataSet_float, _DataSet_float
from pytraj.datasets.DataSet_integer cimport DataSet_integer, _DataSet_integer
from pytraj.datasets.DataSet_string cimport DataSet_string, _DataSet_string
from pytraj.datasets.DataSet_MatrixDbl cimport DataSet_MatrixDbl, _DataSet_MatrixDbl
from pytraj.datasets.DataSet cimport DataSet, _DataSet

def cast_dataset(dsetin, dtype='general'):
    """create memoryview for DataSet instance. 
    DataSet instace is taken from DatatSetList
    Parameters
    ---------
    dset : DataSet instance
    dtype : str (default dtype=None)
        {'general', 'matrix', '1D', '2D', 'double', 'matrix_dbl',
         'integer'}
    """
    # TODO : rename data set
    cdef DataSet dset
    cdef DataSet_1D newset1D
    cdef DataSet_2D newset2D
    cdef DataSet_double newset_double
    cdef DataSet_float newset_float
    cdef DataSet_integer newset_integer
    cdef DataSet_MatrixDbl newset_MatrixDbl
    cdef DataSet_string newset_string

    if not isinstance(dsetin, DataSet):
        dset = <DataSet> dsetin.alloc()
    else:
        dset = <DataSet> dsetin

    dtype = dtype.upper()

    if dtype == '1D':
        newset1D = DataSet_1D()
        # need to recast baseptr0
        newset1D.baseptr0 =  dset.baseptr0
        # need to recast baseptr_1
        newset1D._recast_pointers(0)
        return newset1D

    elif dtype == '2D':
        newset2D = DataSet_2D()
        newset2D.baseptr0 = <_DataSet*> dset.baseptr0
        newset2D.baseptr_1 = <_DataSet_2D*> dset.baseptr0
        return newset2D

    elif dtype in ['GENERAL', 'DOUBLE']: 
        newset_double = DataSet_double()
        # since we introduce memory view, we let cpptraj free memory
        newset_double.py_free_mem = False
        newset_double.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_double.baseptr_1 = <_DataSet_1D*> dset.baseptr0
        newset_double.thisptr = <_DataSet_double*> dset.baseptr0
        return newset_double

    elif dtype in ['FLOAT']: 
        newset_float = DataSet_float()
        # since we introduce memory view, we let cpptraj free memory
        newset_float.py_free_mem = False
        newset_float.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_float.baseptr_1 = <_DataSet_1D*> dset.baseptr0
        newset_float.thisptr = <_DataSet_float*> dset.baseptr0
        return newset_float

    elif dtype in ['INTEGER']: 
        newset_integer = DataSet_integer()
        # since we introduce memory view, we let cpptraj free memory
        newset_integer.py_free_mem = False
        newset_integer.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_integer.baseptr_1 = <_DataSet_1D*> dset.baseptr0
        newset_integer.thisptr = <_DataSet_integer*> dset.baseptr0
        return newset_integer

    elif dtype in ['STRING']: 
        newset_string = DataSet_string()
        # since we introduce memory view, we let cpptraj free memory
        newset_string.py_free_mem = False
        newset_string.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_string.baseptr_1 = <_DataSet_1D*> dset.baseptr0
        newset_string.thisptr = <_DataSet_string*> dset.baseptr0
        return newset_string

    elif dtype in ['MATRIX', 'MATRIX_DBL']:
        newset_matrixdbl = DataSet_MatrixDbl()
        # since we introduce memory view, we let cpptraj free memory
        newset_matrixdbl.py_free_mem = False
        newset_matrixdbl.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_matrixdbl.baseptr_1 = <_DataSet_2D*> dset.baseptr0
        newset_matrixdbl.thisptr = <_DataSet_MatrixDbl*> dset.baseptr0
        return newset_matrixdbl

    else:
        raise NotImplementedError("")
