from __future__ import print_function, absolute_import
import numpy as np
from .c_analysis import c_analysis
from .datasets import c_datasets
from .datasets.c_datasetlist import DatasetList as CpptrajDatasetList

mat_keys = {
    'dist',
    'idea',
    'correl',
    'covar',
    'mwcovar',
    'distcovar',
    'dihcovar',
}

__all__ = mat_keys

__cpptrajdoc__ = """
    cpptraj manual
    --------------
    Calculate a matrix of the specified type from input coordinates.

    + dist: Distance matrix (default).
    + correl: Correlation matrix (aka dynamic cross correlation).
    + covar: Coordinate covariance matrix.
    + mwcovar: Mass-weighted coordinate covariance matrix.
    + distcovar: Distance covariance matrix.
    + idea: Isotropically Distributed Ensemble Analysis matrix.
    + dihcovar: Dihedral covariance matrix.
"""

template = '''
from .c_action import c_action
from .get_common_objects import super_dispatch, get_data_from_dtype

@super_dispatch()
def %s(traj=None, mask="", top=None, dtype='ndarray', mat_type='full', frame_indices=None):
    """Compute matrix

    Parameters
    ----------
    traj : Trajectory-like
    mask : cpptraj mask
    top : Topology, optional, default None
    mat_type : str, {'full', 'half', 'cpptraj'}, default 'full'
        if 'full': 2D full matrix
        if 'half': triangular matrix
        if 'cpptraj': 1D array

    cpptraj compat mode
    -------------------
    {'distance_matrix' : 'dist',
    'correlation_matrix' : 'correl',
    'coord_covariance_matrix' : 'covar',
    'mw_covariance_matrix' : 'mwcovar',
    'distcovar_matrix' : 'distcovar',
    'idea_matrix' : 'idea'}
    """

    dslist = CpptrajDatasetList()
    template_mask = '%s '
    template_mask += mask

    act = c_action.Action_Matrix()
    act(template_mask, traj, top=top, dslist=dslist)
    # need to call `post_process` so cpptraj can normalize some data
    # check cpptraj's code
    act.post_process()
    if dtype == 'ndarray':
        if mat_type == 'full':
            return dslist[0].values
        elif mat_type == 'half':
            return dslist[0].to_half_matrix()
        elif mat_type == 'cpptraj':
            return dslist[0]._to_cpptraj_sparse_matrix()
        else:
            raise ValueError()
    else:
        return get_data_from_dtype(dslist, dtype=dtype)
'''

for k in mat_keys:
    my_func_str = template % (k, k)
    g_dict = globals()
    exec(my_func_str)
    g_dict[k].__doc__ += __cpptrajdoc__

del k

exec('''
from .decorators import register_pmap, register_openmp
dist = register_pmap(dist)
idea = register_pmap(idea)
covar = register_openmp(covar)
''')


def diagonalize(mat, n_vecs, dtype='tuple'):
    '''diagonalize matrix and return (eigenvalues, eigenvectors)

    Parameters
    ----------
    mat : 2D ndarray or DatasetMatrixDouble
    n_vecs : number of output vectors
    dtype : str, {'tuple', 'dataset'}, default 'tuple'
        if 'tuple', return a tuple (eigenvalues, eigenvectors). If 'dataset' return CpptrajDataseList

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> mat = pt.matrix.dist(traj, '@CA')
    >>> x = pt.matrix.diagonalize(mat, 4, dtype='tuple')
    >>> x = pt.matrix.diagonalize(mat, 4, dtype='dataset')

    >>> # use cpptraj dataset to save memory
    >>> mat_cpp = pt.matrix.covar(traj, '@CA', dtype='cpptraj_dataset')[0]
    >>> x = pt.matrix.diagonalize(mat_cpp, 4, dtype='tuple')
    >>> print(x[0].shape, x[1].shape)
    (4,) (4, 36)
    '''
    _vecs = 'vecs ' + str(n_vecs)
    dslist = CpptrajDatasetList()
    dslist.add('matrix_dbl', 'mymat')

    if isinstance(mat, np.ndarray):
        indices = np.triu_indices(mat.shape[0])
        arr = mat[indices]
        dslist[0]._set_data_half_matrix(
            arr.astype('f8'),
            vsize=len(arr),
            n_cols=mat.shape[0])
    elif isinstance(mat, c_datasets.DatasetMatrixDouble):
        dslist[0]._set_data_half_matrix(mat._to_cpptraj_sparse_matrix(),
                                        vsize=mat.size,
                                        n_cols=mat.n_cols)

    act = c_analysis.Analysis_Matrix()
    act(' '.join(('mymat', _vecs)), dslist=dslist)
    dslist._pop(0)

    if dtype == 'tuple':
        return (dslist[-1].eigenvalues, dslist[-1].eigenvectors)
    elif dtype == 'dataset':
        return dslist
    else:
        raise ValueError('only support dtype of tuple or dataset')


def _diagonalize_np(mat, n_vecs):
    '''have opposite sign with cpptraj

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> mat = pt.matrix.dist(traj, '@CA')
    >>> x = pt.matrix._diagonalize_np(mat, 4)
    '''
    evals, evecs = np.linalg.eigh(mat)
    evals = evals[::-1][:n_vecs]
    evecs = evecs[:, ::-1].T[:n_vecs]
    return evals, evecs
