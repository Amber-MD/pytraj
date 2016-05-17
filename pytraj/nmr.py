import numpy as np
from .decorators import register_pmap
from .datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from .c_action import c_action
from .c_analysis import c_analysis
from .c_action.actionlist import ActionList

from .get_common_objects import get_topology, get_data_from_dtype, get_list_of_commands
from .get_common_objects import get_fiterator


def _2darray_to_atommask_groups(seq):
    '''[[0, 3], [4, 7]] turns to ['@1 @4', '@5 @8']
    '''
    for arr in seq:
        # example: arr = [0, 3]; turns ot '@1 @4'
        yield '@' + str(arr[0] + 1) + ' @' + str(arr[1] + 1)


def _ired(iredvec,
          modes,
          NHbond=True,
          freq=None,
          NHdist=1.02,
          order=2,
          tstep=1.0,
          tcorr=10000.,
          norm=False,
          drct=False,
          dtype='dataset'):
    '''perform isotropic reorientational Eigenmode dynamics analysis

    Parameters
    ----------
    iredvec : shape=(n_vectors, n_frames, 3)
    modes : tuple of (eigenvalues, eigenvectors) or DatasetModes (has attribute 'eigenvalues', 'eigenvectors')
        eigenvalues has shape of (n_modes, )
        eigenvectors has shape of (n_modes, vector_size), each vector correspond to each eigenvalue
    NHbond : bool, default True
        if True, freq value will be used
        if False, NHdist and freq will be ignored.
    freq : float or None, default None
        be used with NHbond
    NHdist : N-H bond length, default 1.02
    order : int, default 2
    tstep : timestep between frames, default 1.0 ps
    tcorr: default 10000.
    norm : default False

    Returns
    -------
    CpptrajDatasetList
    '''

    _freq = 'relax freq ' + str(freq) if freq is not None else ''
    _NHdist = 'NHdist ' + str(NHdist) if NHbond else ''
    _order = 'order ' + str(order)
    _tstep = 'tstep ' + str(tstep)
    _tcorr = 'tcorr ' + str(tcorr)
    _norm = 'norm' if norm else ''
    _drct = 'drct' if drct else ''
    _modes = 'modes mymodes'
    if NHbond:
        command = ' '.join((_freq, _NHdist, _modes, _order, _tstep, _tcorr,
                            _norm, _drct))
    else:
        command = ' '.join((_modes, _order, _tstep, _tcorr, _norm, _drct))
    act = c_analysis.Analysis_IRED()
    dslist = CpptrajDatasetList()

    for idx, dvec in enumerate(iredvec):
        name = 'ired_' + str(idx)
        dslist.add('vector', name)
        dslist[-1].scalar_type = 'iredvec'
        dslist[-1].data = np.asarray(dvec, dtype='f8')

    # add data to DatasetModes
    dslist.add('modes', 'mymodes')
    is_reduced = False  # ?
    if hasattr(modes, 'eigenvalues') and hasattr(modes, 'eigenvectors'):
        eigenvalues = modes.eigenvalues
        eigenvectors = modes.eigenvectors
    else:
        eigenvalues, eigenvectors = modes
    dslist[-1]._set_modes(is_reduced, len(eigenvalues), eigenvectors.shape[1],
                          eigenvalues, eigenvectors.flatten())

    act(command, dslist=dslist)
    # remove input datasets to free memory
    # all DatasetVectors + 1 DatasetModes
    # for d in dslist[:idx+2]:
    #    dslist.remove_set(d)
    # values = dslist.values.copy()
    # return values
    # return get_data_from_dtype(dslist, dtype=dtype)
    return dslist


@register_pmap
def ired_vector_and_matrix(traj=None,
                           mask="",
                           frame_indices=None,
                           order=2,
                           dtype='tuple',
                           top=None):
    """perform vector calculation and then calculate ired matrix

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces :class:`pytraj.Frame`
    mask : str or a list of strings
    frame_indices : array-like, optional, default None
        only perform calculation for given frame indices
    order : default 2
    dtype : output's dtype, {'dataset', 'tuple'} default 'dataset'
    top : Topology, optional, default None

    Returns
    -------
    out : if dtype is 'dataset', return pytraj.DatasetList with shape=(n_vectors+1,)
        last index is a matrix, otherwise n_vectors. If dtype is 'tuple', return a a tuple
        of vectors and matrix

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> h = pt.select('@H', traj.top)
    >>> n = h - 1
    >>> nh = list(zip(n, h))
    >>> vecs, mat = pt.ired_vector_and_matrix(traj, mask=nh)
    >>> dslist = pt.ired_vector_and_matrix(traj, mask=nh, dtype='dataset')

    """
    dslist = CpptrajDatasetList()
    _top = get_topology(traj, top)
    fi = get_fiterator(traj, frame_indices)

    cm_arr = np.asarray(mask)

    if cm_arr.dtype.kind != 'i':
        list_of_commands = get_list_of_commands(mask)
    else:
        if cm_arr.ndim != 2:
            raise ValueError(
                'if mask is a numpy.ndarray, it must have ndim = 2')
        list_of_commands = _2darray_to_atommask_groups(cm_arr)

    actlist = ActionList()

    for command in list_of_commands:
        # tag ired vector
        command += ' ired '
        act = c_action.Action_Vector()
        actlist.add(act, command, _top, dslist=dslist)

    act_matired = c_action.Action_Matrix()
    ired_cm = 'ired order ' + str(order)
    actlist.add(act_matired, ired_cm, _top, dslist=dslist)

    actlist.compute(fi)

    if dtype == 'tuple':
        mat = dslist[-1].values
        mat = mat / mat[0, 0]
        return dslist[:-1].values, mat
    else:
        out = get_data_from_dtype(dslist, dtype=dtype)
        if dtype == 'dataset':
            out[-1].values = out[-1].values / out[-1].values[0, 0]
        else:
            out = out
        return out


calc_ired_vector_and_matrix = ired_vector_and_matrix


@register_pmap
def NH_order_parameters(traj,
                        vector_pairs,
                        order=2,
                        tstep=1.,
                        tcorr=10000.,
                        n_cores=1, **kwargs):
    '''compute NH order parameters

    Parameters
    ----------
    traj : Trajectory-like
    vector_pairs : 2D array-like, shape (n_pairs, 2)
    order : default 2
    tstep : default 1.
    tcorr : default 10000.
    kwargs : additional keyword argument
 
    Returns
    -------
    S2 : 1D array, order parameters

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> h_indices = pt.select('@H', traj.top)
    >>> n_indices = h_indices - 1
    >>> nh_pairs = list(zip(n_indices, h_indices))
    >>> data = pt.NH_order_parameters(traj, nh_pairs)
    '''
    from pytraj import matrix

    # compute N-H vectors and ired matrix
    if n_cores == 1:
        vecs_and_mat = ired_vector_and_matrix(traj,
                                              vector_pairs,
                                              order=order,
                                              dtype='tuple')
    else:
        # use _pmap to avoid cicular import
        from pytraj import _pmap
        vecs_and_mat = _pmap(ired_vector_and_matrix,
                             traj,
                             vector_pairs,
                             order=order,
                             dtype='tuple',
                             n_cores=n_cores)

    state_vecs = vecs_and_mat[0]
    mat_ired = vecs_and_mat[1]

    # get eigenvalues and eigenvectors
    modes = matrix.diagonalize(mat_ired, n_vecs=len(state_vecs), dtype='dataset')[0]
    evals, evecs = modes.eigenvalues, modes.eigenvectors

    data = _ired(state_vecs,
                 modes=(evals, evecs),
                 NHbond=True,
                 tcorr=tcorr,
                 tstep=tstep, **kwargs)
    order = [d.values.copy() for d in data if 'S2' in d.key][0]
    return order

# create alis for easy parsing
calc_NH_order_parameters = NH_order_parameters
