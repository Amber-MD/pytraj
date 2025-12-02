"""
Principal Component Analysis functions
"""
from .base import *

__all__ = ['pca', 'projection']


@super_dispatch()
def projection(traj,
               mask='',
               eigenvectors=None,
               eigenvalues=None,
               scalar_type='covar',
               average_coords=None,
               frame_indices=None,
               dtype='ndarray',
               top=None):
    '''compute projection along given eigenvectors

    Parameters
    ----------
    traj : Trajectory-like
    mask : atom mask, either string or array-like
    eigenvalues : 1D array-like
    eigenvectors : 2D array-like
    scalar_type : str, {'covar', 'mwcovar', }, default 'covar'
        make sure to provide correct scalar_type.
        Note: not yet support 'dihcovar' and 'idea'
    average_coords : 3D array-like, optional, default None
        average coordinates. If 'None', pytraj will compute mean_structure with given mask
    frame_indices : array-like
        If not None, compute projection for given frames.
    dtype : str, return data type, default 'ndarray'
    top : Topology, optional, default None

    Returns
    -------
    projection_data : ndarray, shape=(n_vecs, n_frames)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> mat = pt.matrix.covar(traj, '@CA')
    >>> eigenvectors, eigenvalues = pt.matrix.diagonalize(mat, 2)

    >>> # since we compute covariance matrix, we need to specify
    >>> # scalar_type = 'covar'
    >>> scalar_type = 'covar'
    >>> data = pt.projection(traj, '@CA', eigenvalues=eigenvalues, eigenvectors=eigenvectors, scalar_type=scalar_type)
    >>> data.shape
    (2, 101)
    '''
    from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList
    from .base import DatasetType
    from ..analysis.c_action import c_action
    from .utilities import mean_structure
    from ..utils.get_common_objects import get_data_from_dtype

    projection_action = c_action.Action_Projection()
    action_datasets = CpptrajDatasetList()

    mode_name = 'my_modes'
    action_datasets.add(DatasetType.MODES, mode_name)

    is_reduced = False
    dataset_mode = action_datasets[-1]
    n_vectors = len(eigenvalues)
    dataset_mode._set_modes(is_reduced, n_vectors, eigenvectors.shape[1],
                            eigenvalues, eigenvectors.flatten())
    dataset_mode.scalar_type = scalar_type

    if average_coords is None:
        frame = mean_structure(traj, mask)
        average_coords = frame.xyz

    dataset_mode._allocate_avgcoords(3 * average_coords.shape[0])
    dataset_mode._set_avg_frame(average_coords.flatten())

    command = f"evecs {mode_name} {mask} beg 1 end {n_vectors}"
    projection_action(command, traj, top=top, dslist=action_datasets)

    action_datasets._pop(0)

    return get_data_from_dtype(action_datasets, dtype=dtype)


def pca(traj,
        mask,
        n_vecs=2,
        fit=True,
        ref=None,
        ref_mask=None,
        dtype='ndarray',
        top=None):
    '''perform PCA analysis by following below steps:

    - (optional) perform rmsfit to reference if needed
    - compute covariance matrix
    - diagonalize the matrix to get eigenvectors and eigenvalues
    - perform projection of each frame with mask to each eigenvector

    Parameters
    ----------
    traj : Trajectory
        traj must be ``pytraj.Trajectory``, which can be created by ``pytraj.load`` method.
    mask : str
        atom mask for covariance matrix and projection
    n_vecs : int, default 2
        number of eigenvectors. If user want to compute projection for all eigenvectors,
        should specify n_vecs=-1 (or a negative number)
    fit : bool, default True
        if True, perform fitting before compute covariance matrix
        if False, no fitting (keep rotation and translation). In this case, `pytraj` will ignore `ref` argument.
    ref : {None, Frame, int}, default None
        if None, trajectory will be superposed to average structure
        if is Frame or integer value, trajectory will be superposed to given reference
    ref_mask : {None, str}, default None (use `mask`)
        if None, use `mask` for fitting
        if str, use this given mask for fitting
    dtype : return datatype
    top : Topology, optional

    Returns
    -------
    out1: projection_data, ndarray with shape=(n_vecs, n_frames)
    out2: tuple of (eigenvalues, eigenvectors)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()[:]

    >>> # compute pca for first and second modes
    >>> pca_data = pt.pca(traj, '!@H=', n_vecs=2)
    >>> # get projection values
    >>> pca_data[0] # doctest: +SKIP
    array([[  4.93425131,  13.80002308,  20.61605835, ..., -57.92280579,
            -61.25728607, -52.85142136],
           [  4.03333616,  -6.9132452 , -14.53991318, ...,  -6.757936  ,
              2.1086719 ,  -3.60922861]], dtype=float32)
    >>> # get eigenvalues for first 2 modes
    >>> pca_data[1][0] # doctest: +SKIP
    array([ 1399.36472919,   240.42342439])

    >>> # compute pca for all modes
    >>> pca_data = pt.pca(traj, '!@H=', n_vecs=-1)

    >>> # does not perform fitting
    >>> data = pt.pca(traj, mask='!@H=', fit=False)

    >>> # provide different mask for fitting
    >>> data = pt.pca(traj, mask='!@H=', fit=True, ref=0, ref_mask='@CA')

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()[:]

    >>> # compute pca for first and second modes
    >>> pca_data = pt.pca(traj, '!@H=', n_vecs=2)
    >>> # get projection values
    >>> pca_data[0] # doctest: +SKIP
    array([[  4.93425131,  13.80002308,  20.61605835, ..., -57.92280579,
            -61.25728607, -52.85142136],
           [  4.03333616,  -6.9132452 , -14.53991318, ...,  -6.757936  ,
              2.1086719 ,  -3.60922861]], dtype=float32)
    >>> # get eigenvalues for first 2 modes
    >>> pca_data[1][0] # doctest: +SKIP
    array([ 1399.36472919,   240.42342439])

    >>> # compute pca for all modes
    >>> pca_data = pt.pca(traj, '!@H=', n_vecs=-1)

    >>> # does not perform fitting
    >>> data = pt.pca(traj, mask='!@H=', fit=False)

    >>> # provide different mask for fitting
    >>> data = pt.pca(traj, mask='!@H=', fit=True, ref=0, ref_mask='@CA')
    '''
    from ..trajectory.trajectory import Trajectory
    from ..trajectory.trajectory_iterator import TrajectoryIterator
    from ..analysis import matrix
    from .utilities import mean_structure
    from ..utils.get_common_objects import get_reference

    ref_mask = ref_mask if ref_mask is not None else mask

    if not isinstance(traj, (Trajectory, TrajectoryIterator)):
        raise ValueError('must be Trajectory-like')

    if fit:
        if ref is None:
            traj.superpose(ref=0, mask=ref_mask)
            avg = mean_structure(traj)
            traj.superpose(ref=avg, mask=ref_mask)
            n_refs = 2
        else:
            ref = get_reference(traj, ref)
            traj.superpose(ref=ref, mask=ref_mask)
            n_refs = 1

    avg2 = mean_structure(traj, mask=mask)

    covariance_matrix = matrix.covar(traj, mask)
    n_vecs = covariance_matrix.shape[0] if n_vecs < 0 else n_vecs

    eigenvectors, eigenvalues = matrix.diagonalize(covariance_matrix,
                                                   n_vecs=n_vecs,
                                                   dtype='tuple')

    projection_data = projection(traj,
                                 mask=mask,
                                 average_coords=avg2.xyz,
                                 eigenvalues=eigenvalues,
                                 eigenvectors=eigenvectors,
                                 scalar_type='covar',
                                 dtype=dtype)

    if fit and hasattr(traj, '_transform_commands'):
        for _ in range(n_refs):
            traj._transform_commands.pop()
        if traj._transform_commands:
            traj._reset_transformation()
        else:
            traj._remove_transformations()

    return projection_data, (eigenvalues, eigenvectors)
