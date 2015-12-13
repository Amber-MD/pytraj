# distutils: language = c++
from .c_datasets cimport(_Dataset, Dataset, Dataset1D, _Dataset1D, DatasetInteger, _DatasetInteger,
                           _DatasetFloat, DatasetFloat, DatasetDouble, _DatasetDouble,
                           DatasetString, _DatasetString, _DatasetVector, DatasetVector,
                           Dataset2D, _Dataset2D, DatasetMatrixDouble,
                           _DatasetMatrixDouble, DatasetMatrixFloat, _DatasetMatrixFloat,
                           Dataset3D, _Dataset3D, DatasetGridFloat, _DatasetGridFloat,
                           _DatasetModes, DatasetModes,
                           DatasetMesh, _DatasetMesh, _DatasetMatrix3x3, DatasetMatrix3x3,
                           _DatasetCoords, DatasetCoords, _DatasetCoordsRef,
                           DatasetCoordsRef, _DatasetCoordsCRD, DatasetCoordsCRD,
                           DatasetTopology, _DatasetTopology)
from ..c_traj.c_trajectory cimport TrajectoryCpptraj, _TrajectoryCpptraj


def cast_dataset(dsetin=None, dtype='general'):
    """cast cpptraj's Dataset to other
    Dataset instace is taken from DatatSetList

    Parameters
    ---------
    dset : Dataset instance
    dtype : str (default dtype=None)
        {'general', 'matrix', '1D', '2D', 'double',
         'mesh',
         'matrix_dbl', 'matrix_flt',
         'integer',
         'coords_crd', 'coords'
         'coords_trj', 'trj',
         'topology',}
    """
    # TODO:
    cdef Dataset dset
    cdef Dataset1D newset1D
    cdef Dataset2D newset2D
    cdef DatasetDouble newset_double
    cdef DatasetFloat newset_float
    cdef DatasetInteger newset_integer
    cdef DatasetString newset_string
    cdef DatasetModes newset_modes
    cdef DatasetMesh newset_mesh
    cdef DatasetVector newset_vector
    cdef DatasetMatrix3x3 newset_matrix3x3
    cdef DatasetMatrixDouble newset_matrixdbl
    cdef DatasetMatrixFloat newset_matrixflt
    cdef DatasetGridFloat newset_gridflt
    cdef DatasetCoordsRef newset_coords_ref
    cdef DatasetCoordsCRD newset_coords_crd
    cdef DatasetTopology newset_topology
#    cdef TrajectoryCpptraj newset_coords_trj

    if not isinstance(dsetin, Dataset):
        dset = <Dataset> dsetin.alloc()
    else:
        dset = <Dataset> dsetin

    dtype = dtype.upper()

    if dtype == '1D':
        newset1D = Dataset1D()
        # need to recast baseptr0
        newset1D.baseptr0 = dset.baseptr0
        # need to recast baseptr_1
        newset1D._recast_pointers(0)
        return newset1D

    elif dtype == '2D':
        newset2D = Dataset2D()
        newset2D.baseptr0 = <_Dataset*> dset.baseptr0
        newset2D.baseptr_1 = <_Dataset2D*> dset.baseptr0
        return newset2D

    elif dtype in ['GENERAL', 'DOUBLE']:
        newset_double = DatasetDouble()
        # since we introduce memory view, we let cpptraj free memory
        newset_double._own_memory = False
        newset_double.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_double.baseptr_1 = <_Dataset1D*> dset.baseptr0
        newset_double.thisptr = <_DatasetDouble*> dset.baseptr0
        return newset_double

    elif dtype in ['FLOAT']:
        newset_float = DatasetFloat()
        # since we introduce memory view, we let cpptraj free memory
        newset_float._own_memory = False
        newset_float.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_float.baseptr_1 = <_Dataset1D*> dset.baseptr0
        newset_float.thisptr = <_DatasetFloat*> dset.baseptr0
        return newset_float

    elif dtype in ['INTEGER']:
        newset_integer = DatasetInteger()
        # since we introduce memory view, we let cpptraj free memory
        newset_integer._own_memory = False
        newset_integer.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_integer.baseptr_1 = <_Dataset1D*> dset.baseptr0
        newset_integer.thisptr = <_DatasetInteger*> dset.baseptr0
        return newset_integer

    elif dtype in ['STRING']:
        newset_string = DatasetString()
        # since we introduce memory view, we let cpptraj free memory
        newset_string._own_memory = False
        newset_string.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_string.baseptr_1 = <_Dataset1D*> dset.baseptr0
        newset_string.thisptr = <_DatasetString*> dset.baseptr0
        return newset_string

    elif dtype in ['XYMESH', 'MESH']:
        newset_mesh = DatasetMesh()
        # since we introduce memory view, we let cpptraj free memory
        newset_mesh._own_memory = False
        newset_mesh.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_mesh.baseptr_1 = <_Dataset1D*> dset.baseptr0
        newset_mesh.thisptr = <_DatasetMesh*> dset.baseptr0
        return newset_mesh

    elif dtype in ['MODES']:
        newset_modes = DatasetModes()
        # since we introduce memory view, we let cpptraj free memory
        newset_modes._own_memory = False
        newset_modes.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_modes.thisptr = <_DatasetModes*> dset.baseptr0
        return newset_modes

    elif dtype in ['VECTOR']:
        newset_vector = DatasetVector()
        # since we introduce memory view, we let cpptraj free memory
        newset_vector._own_memory = False
        newset_vector.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_vector.thisptr = <_DatasetVector*> dset.baseptr0
        return newset_vector

    elif dtype in ['MAT3X3', 'MATRIX3X3']:
        newset_matrix3x3 = DatasetMatrix3x3()
        # since we introduce memory view, we let cpptraj free memory
        newset_matrix3x3._own_memory = False
        newset_matrix3x3.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_matrix3x3.thisptr = <_DatasetMatrix3x3*> dset.baseptr0
        return newset_matrix3x3

    elif dtype in ['MATRIX_DBL', 'MATRIX_DOUBLE', 'MATRIX DOUBLE']:
        newset_matrixdbl = DatasetMatrixDouble()
        # since we introduce memory view, we let cpptraj free memory
        newset_matrixdbl._own_memory = False
        newset_matrixdbl.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_matrixdbl.baseptr_1 = <_Dataset2D*> dset.baseptr0
        newset_matrixdbl.thisptr = <_DatasetMatrixDouble*> dset.baseptr0
        return newset_matrixdbl

    elif dtype in ['MATRIX_FLT', 'MATRIX_FLOAT', 'MATRIX FLOAT']:
        newset_matrixflt = DatasetMatrixFloat()
        # since we introduce memory view, we let cpptraj free memory
        newset_matrixflt._own_memory = False
        newset_matrixflt.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_matrixflt.baseptr_1 = <_Dataset2D*> dset.baseptr0
        newset_matrixflt.thisptr = <_DatasetMatrixFloat*> dset.baseptr0
        return newset_matrixflt

    elif dtype in ['GRID_FLT', 'GRID_FLOAT', 'GRID FLOAT', 'GRID']:
        newset_gridflt = DatasetGridFloat()
        # since we introduce memory view, we let cpptraj free memory
        newset_gridflt._own_memory = False
        newset_gridflt.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_gridflt.baseptr_1 = <_Dataset3D*> dset.baseptr0
        newset_gridflt.thisptr = <_DatasetGridFloat*> dset.baseptr0
        return newset_gridflt

    elif dtype in ['COORDS_CRD', 'COORDS', 'CRD']:
        # FIXME: not correctly casting
        # get '0' size when casting back from Dataset
        newset_coords_crd = DatasetCoordsCRD()
        # since we introduce memory view, we let cpptraj free memory
        newset_coords_crd._own_memory = False
        newset_coords_crd.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_coords_crd.baseptr_1 = <_DatasetCoords*> dset.baseptr0
        newset_coords_crd.thisptr = <_DatasetCoordsCRD*> dset.baseptr0
        return newset_coords_crd

    elif dtype in ['COORDS_TRJ', 'TRJ', 'TRAJ', 'COORDS_TRAJ']:
        newset_coords_trj = TrajectoryCpptraj()
        # since we introduce memory view, we let cpptraj free memory
        newset_coords_trj._own_memory = False
        # make sure other pointers pointing to the same address
        newset_coords_trj.thisptr = <_TrajectoryCpptraj*> dset.baseptr0
        return newset_coords_trj

    elif dtype in ['COORDS_REF_FRAME', 'REF_FRAME', 'REFFRAME', 'REF', 'REFERENCE']:
        newset_coords_ref = DatasetCoordsRef()
        # since we introduce memory view, we let cpptraj free memory
        newset_coords_ref._own_memory = False
        newset_coords_ref.baseptr0 = dset.baseptr0
        # make sure other pointers pointing to the same address
        newset_coords_ref.baseptr_1 = <_DatasetCoords*> dset.baseptr0
        newset_coords_ref.thisptr = <_DatasetCoordsRef*> dset.baseptr0
        return newset_coords_ref
    elif dtype in ['TOPOLOGY']:
        newset_topology = DatasetTopology()
        newset_topology._own_memory = False
        newset_topology.baseptr0 = dset.baseptr0
        newset_topology.thisptr = <_DatasetTopology*> dset.baseptr0
        return newset_topology
    else:
        raise NotImplementedError("")
