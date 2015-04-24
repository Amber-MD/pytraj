# distutils: language = c++
from __future__ import absolute_import

from .DataSet_1D cimport DataSet_1D, _DataSet_1D
from .DataSet_2D cimport DataSet_2D, _DataSet_2D
from .DataSet_3D cimport DataSet_3D, _DataSet_3D
from .DataSet_double cimport DataSet_double, _DataSet_double
from .DataSet_float cimport DataSet_float, _DataSet_float
from .DataSet_integer cimport DataSet_integer, _DataSet_integer
from .DataSet_string cimport DataSet_string, _DataSet_string
from .DataSet_Mesh cimport DataSet_Mesh, _DataSet_Mesh
from .DataSet_Vector cimport _DataSet_Vector, DataSet_Vector
from .DataSet_MatrixDbl cimport DataSet_MatrixDbl, _DataSet_MatrixDbl
from .DataSet_MatrixFlt cimport DataSet_MatrixFlt, _DataSet_MatrixFlt
from .DataSet_GridFlt cimport DataSet_GridFlt, _DataSet_GridFlt
from .DataSet cimport DataSet, _DataSet
from .DataSet_Coords cimport _DataSet_Coords, DataSet_Coords
from .DataSet_Coords_REF cimport _DataSet_Coords_REF, DataSet_Coords_REF
from .DataSet_Coords_CRD cimport _DataSet_Coords_CRD, DataSet_Coords_CRD
from .DataSet_Coords_TRJ cimport _DataSet_Coords_TRJ, DataSet_Coords_TRJ

def cast_dataset(dsetin=None, dtype='general'):
    """create memoryview for DataSet instance. 
    DataSet instace is taken from DatatSetList
    Parameters
    ---------
    dset : DataSet instance
    dtype : str (default dtype=None)
        {'general', 'matrix', '1D', '2D', 'double', 
         'mesh',
         'matrix_dbl', 'matrix_flt',
         'integer',
         'coords_crd', 'coords'
         'coords_trj', 'trj'}
    """
    # TODO:
    cdef DataSet dset
    cdef DataSet_1D newset1D
    cdef DataSet_2D newset2D
    cdef DataSet_double newset_double
    cdef DataSet_float newset_float
    cdef DataSet_integer newset_integer
    cdef DataSet_string newset_string
    cdef DataSet_Mesh newset_mesh
    cdef DataSet_Vector newset_vector
    cdef DataSet_MatrixDbl newset_matrixdbl
    cdef DataSet_MatrixFlt newset_matrixflt
    cdef DataSet_GridFlt newset_gridflt
    cdef DataSet_Coords_REF newset_coords_ref
    cdef DataSet_Coords_CRD newset_coords_crd
    cdef DataSet_Coords_TRJ newset_coords_trj

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

    elif dtype in ['XYMESH']: 
        newset_mesh = DataSet_Mesh()
        # since we introduce memory view, we let cpptraj free memory
        newset_mesh.py_free_mem = False
        newset_mesh.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_mesh.baseptr_1 = <_DataSet_1D*> dset.baseptr0
        newset_mesh.thisptr = <_DataSet_Mesh*> dset.baseptr0
        return newset_mesh

    elif dtype in ['VECTOR']: 
        newset_vector = DataSet_Vector()
        # since we introduce memory view, we let cpptraj free memory
        newset_vector.py_free_mem = False
        newset_vector.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_vector.baseptr_1 = <_DataSet_1D*> dset.baseptr0
        newset_vector.thisptr = <_DataSet_Vector*> dset.baseptr0
        return newset_vector

    elif dtype in ['MATRIX_DBL', 'MATRIX_DOUBLE', 'MATRIX DOUBLE']:
        newset_matrixdbl = DataSet_MatrixDbl()
        # since we introduce memory view, we let cpptraj free memory
        newset_matrixdbl.py_free_mem = False
        newset_matrixdbl.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_matrixdbl.baseptr_1 = <_DataSet_2D*> dset.baseptr0
        newset_matrixdbl.thisptr = <_DataSet_MatrixDbl*> dset.baseptr0
        return newset_matrixdbl

    elif dtype in ['MATRIX_FLT', 'MATRIX_FLOAT', 'MATRIX FLOAT']:
        newset_matrixflt = DataSet_MatrixFlt()
        # since we introduce memory view, we let cpptraj free memory
        newset_matrixflt.py_free_mem = False
        newset_matrixflt.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_matrixflt.baseptr_1 = <_DataSet_2D*> dset.baseptr0
        newset_matrixflt.thisptr = <_DataSet_MatrixFlt*> dset.baseptr0
        return newset_matrixflt

    elif dtype in ['GRID_FLT', 'GRID_FLOAT', 'GRID FLOAT']:
        newset_gridflt = DataSet_GridFlt()
        # since we introduce memory view, we let cpptraj free memory
        newset_gridflt.py_free_mem = False
        newset_gridflt.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_gridflt.baseptr_1 = <_DataSet_3D*> dset.baseptr0
        newset_gridflt.thisptr = <_DataSet_GridFlt*> dset.baseptr0
        return newset_gridflt

    elif dtype in ['COORDS_CRD', 'COORDS', 'CRD']:
        # FIXME: not correctly casting
        # get '0' size when casting back from DataSet
        newset_coords_crd = DataSet_Coords_CRD()
        # since we introduce memory view, we let cpptraj free memory
        newset_coords_crd.py_free_mem = False
        newset_coords_crd.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_coords_crd.baseptr_1 = <_DataSet_Coords*> dset.baseptr0
        newset_coords_crd.thisptr = <_DataSet_Coords_CRD*> dset.baseptr0
        return newset_coords_crd

    elif dtype in ['COORDS_TRJ', 'TRJ', 'TRAJ', 'COORDS_TRAJ']:
        newset_coords_trj = DataSet_Coords_TRJ()
        # since we introduce memory view, we let cpptraj free memory
        newset_coords_trj.py_free_mem = False
        newset_coords_trj.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_coords_trj.baseptr_1 = <_DataSet_Coords*> dset.baseptr0
        newset_coords_trj.thisptr = <_DataSet_Coords_TRJ*> dset.baseptr0
        return newset_coords_trj

    elif dtype in ['COORDS_REF_FRAME', 'REF_FRAME', 'REFFRAME', 'REF', 'REFERENCE']:
        newset_coords_ref = DataSet_Coords_REF()
        # since we introduce memory view, we let cpptraj free memory
        newset_coords_ref.py_free_mem = False
        newset_coords_ref.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_coords_ref.baseptr_1 = <_DataSet_Coords*> dset.baseptr0
        newset_coords_ref.thisptr = <_DataSet_Coords_REF*> dset.baseptr0
        return newset_coords_ref
    else:
        raise NotImplementedError("")
