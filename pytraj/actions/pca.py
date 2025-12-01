"""
Principal Component Analysis functions
"""
from .base import *

__all__ = [
    'pca', 'projection'
]


@super_dispatch()
def projection(traj,
               mask='',
               mode=0,
               beg=1,
               end=-1,
               eigenfile='',
               efile='',
               dtype='ndarray',
               scalar_type='covar',
               options='',
               top=None,
               frame_indices=None):
    """compute projection

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    mode : int, default 0
        mode index
    beg : int, default 1
        first mode
    end : int, default -1
        last mode
    eigenfile : str, optional
        eigenvalue file
    efile : str, optional
        eigenvector file
    dtype : str, default 'ndarray'
        return data type
    scalar_type : str, default 'covar'
        scalar type (covar, mwcovar, correl, distcovar, idea, ired)
    options : str, optional
        extra options
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray or DatasetList
    """
    command = mask
    if mode != 0:
        command += f" mode {mode}"
    if beg != 1:
        command += f" beg {beg}"
    if end != -1:
        command += f" end {end}"
    if eigenfile:
        command += f" eigenfile {eigenfile}"
    if efile:
        command += f" efile {efile}"
    command += f" {scalar_type} {options}"

    c_dslist = CpptrajDatasetList()
    c_action = c_action.Action_Projection()
    c_action.read_input(command, top=traj.top, dslist=c_dslist)
    c_action.setup(traj.top)

    for frame in traj:
        c_action.compute(frame)

    c_action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


def pca(traj,
        mask='*',
        n_vecs=2,
        dtype='ndarray',
        top=None,
        frame_indices=None):
    """perform principle component analysis

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*' (all atoms)
        atom mask
    n_vecs : int, default 2
        number of eigenvectors to save
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : DatasetList if dtype is 'dataset', otherwise return ndarray

    Notes
    -----
    return eigenvalues and eigenvectors

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.pca(traj, mask='!@H=', n_vecs=2)
    >>> data[0].shape
    (2,)
    >>> data[1].shape
    (2, 7998)
    """
    # use matrix for PCA calculation
    from .matrix import matrix

    # get covariance matrix
    mat = matrix(traj, mask=mask, byatom=True, dtype='dataset')
    covar_mat = mat['matrix_covar']

    # perform eigenvalue decomposition
    c_dslist = CpptrajDatasetList()

    # add covariance matrix to dataset list
    mat_dataset = c_dslist.add('matrix_dbl', 'covar')
    mat_dataset.data = covar_mat

    command = f"covar out evecs.dat vecs {n_vecs}"

    # run analysis
    c_analysis.Analysis_Matrix(command, dslist=c_dslist)

    if dtype == 'dataset':
        return c_dslist
    else:
        # extract eigenvalues and eigenvectors
        eigenvals = None
        eigenvecs = None

        for dataset in c_dslist:
            if 'eigenvalues' in dataset.legend:
                eigenvals = dataset.values
            elif 'eigenvectors' in dataset.legend or 'evecs' in dataset.legend:
                eigenvecs = dataset.values

        if eigenvals is not None and eigenvecs is not None:
            return [eigenvals, eigenvecs]
        elif eigenvals is not None:
            return eigenvals
        elif eigenvecs is not None:
            return eigenvecs
        else:
            return get_data_from_dtype(c_dslist, dtype=dtype)