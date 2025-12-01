"""
Geometric analysis functions: distances, angles, dihedrals
"""
from .base import *

__all__ = [
    'distance', 'pairwise_distance', 'angle', 'dihedral', 'mindist',
    'dihedral_rms', 'distance_to_point', 'distance_to_reference', '_calculate_distance',
    '_calculate_angles_for_int_array', '_calculate_dihedrals_for_int_array', '_dihedral_res'
]


def _calculate_distance(traj, int_2darr: np.ndarray, n_frames: int, dtype: str) -> Union[np.ndarray, DatasetList]:
    if int_2darr.ndim == 1:
        int_2darr = np.atleast_2d(int_2darr)

    arr = np.empty([n_frames, len(int_2darr)])

    for idx, frame in enumerate(iterframe_master(traj)):
        arr[idx] = frame._distance(int_2darr)

    arr = arr.T
    if dtype == 'ndarray':
        return arr
    else:
        dslist = DatasetList({'distance': arr})
        return get_data_from_dtype(dslist, dtype)


def _calculate_angles_for_int_array(traj, integer_array, n_frames, dtype):
    if integer_array.ndim == 1:
        integer_array = np.atleast_2d(integer_array)

    if integer_array.shape[1] != 3:
        raise ValueError("require int-array with shape=(n_atoms, 3)")

    arr = np.empty([n_frames, len(integer_array)])
    for idx, frame in enumerate(iterframe_master(traj)):
        arr[idx] = frame._angle(integer_array)

    arr = arr.T
    if dtype == 'ndarray':
        return arr
    else:
        py_dslist = DatasetList({'angle': arr})
        return get_data_from_dtype(py_dslist, dtype)


def _calculate_dihedrals_for_int_array(traj, integer_array, n_frames, dtype):
    if integer_array.ndim == 1:
        integer_array = np.atleast_2d(integer_array)

    if integer_array.shape[1] != 4:
        raise ValueError("require int-array with shape=(n_atoms, 4)")

    arr = np.empty([n_frames, len(integer_array)])
    for idx, frame in enumerate(iterframe_master(traj)):
        arr[idx] = frame._dihedral(integer_array)

    arr = arr.T
    if dtype == 'ndarray':
        return arr
    else:
        py_dslist = DatasetList({'dihedral': arr})
        return get_data_from_dtype(py_dslist, dtype)


def _dihedral_res(traj, mask=(), resid=0, dtype='ndarray', top=None):
    '''compute dihedral within a single residue. For internal use only.

    Parameters
    ----------
    traj : Trajectory-like
    mask : tuple of strings
    resid : str, resid
    dtype

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.all_actions import _dihedral_res
    >>> traj = pt.datafiles.load_tz2()
    >>> data = _dihedral_res(traj, mask=('N', 'CA', 'C', 'O'), resid=0)
    >>> # use string for resid
    >>> data = _dihedral_res(traj, mask=('N', 'CA', 'C', 'O'), resid='1')
    '''

    if is_int(resid):
        resid = str(resid + 1)
    else:
        resid = resid
    m = ' :%s@' % resid
    command = m + m.join(mask)
    return dihedral(traj=traj, mask=command, top=top, dtype=dtype)


def distance_to_point(traj, mask, point, **kwargs):
    """Calculate distance from atoms to a point

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    point : array-like or tuple
        coordinates of the point
    **kwargs : additional arguments
        passed to _distance_to_ref_or_point

    Returns
    -------
    ndarray : distances from atoms to point
    """
    return _distance_to_ref_or_point(traj=traj, mask=mask, point=point, **kwargs)


def distance_to_reference(traj, mask, ref, **kwargs):
    """Calculate distance from atoms to a reference

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    ref : Trajectory-like or Frame
        reference structure
    **kwargs : additional arguments
        passed to _distance_to_ref_or_point

    Returns
    -------
    ndarray : distances from atoms to reference
    """
    return _distance_to_ref_or_point(traj=traj, mask=mask, ref=ref, **kwargs)


def _calculate_distance(traj, int_2darr: np.ndarray, n_frames: int, dtype: str) -> Union[np.ndarray, DatasetList]:
    """Helper function to calculate distances"""

    if dtype != 'ndarray':
        # Use standard pytraj approach for other dtypes
        dslist = CpptrajDatasetList()

        for idx in int_2darr:
            command = f"distance :{idx[0]+1} :{idx[1]+1}"
            _distance = c_action.Action_Distance()
            _distance.read_input(command, top=traj.top, dslist=dslist)
            _distance.setup(traj.top)

            for frame in traj:
                _distance.compute(frame)

            _distance.post_process()

        return get_data_from_dtype(dslist, dtype=dtype)
    else:
        # Direct numpy calculation for speed
        distances = np.empty((len(int_2darr), n_frames))
        for i, (idx0, idx1) in enumerate(int_2darr):
            traj_xyz_subset = traj.xyz[:, [idx0, idx1]]  # shape: (n_frames, 2, 3)
            diff = traj_xyz_subset[:, 0] - traj_xyz_subset[:, 1]  # shape: (n_frames, 3)
            distances[i] = np.sqrt(np.sum(diff**2, axis=1))  # shape: (n_frames,)

        if len(int_2darr) == 1:
            return distances[0]
        else:
            return distances


def distance(traj=None,
             command='',
             mask='',
             indices=None,
             top=None,
             dtype='ndarray',
             frame_indices=None):
    """Compute distance.

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        cpptraj distance command
    mask : {str, array-like}, optional
        Atom mask(s)
    indices : {None, array-like}, optional
        If given, use array indices. If not given, use `mask`
    top : Topology, optional, default=None
    dtype : str, default='ndarray'
        return data type
    frame_indices : array-like, optional
        frame indices

    Returns
    -------
    distances : ndarray or list of 1D ndarray if has more than one distance

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # single distance, using masks
    >>> data = pt.distance(traj, mask=['@CA', '@CB'])
    >>> data.shape
    (10,)

    >>> # multiple distances, using masks
    >>> data = pt.distance(traj, mask=[['@CA', '@CB'], ['@N', '@H']])
    >>> data.shape
    (2, 10)

    >>> # using cpptraj command
    >>> data = pt.distance(traj, command="distance :1 :2")
    >>> data.shape
    (10,)
    """

    top = traj.top.copy() if top is None else top

    if indices is not None:
        int_2darr = np.asarray(indices, dtype=int)
        if int_2darr.ndim != 2:
            raise ValueError("indices must be 2D array")
        return _calculate_distance(traj, int_2darr, traj.n_frames, dtype)

    elif mask:
        indices = get_list_of_commands(mask, top=top)
        int_2darr = np.asarray(indices, dtype=int)
        return _calculate_distance(traj, int_2darr, traj.n_frames, dtype)

    elif command:
        dslist = CpptrajDatasetList()
        act = c_action.Action_Distance()
        act.read_input(command, top=top, dslist=dslist)
        act.setup(top)

        for frame in traj:
            act.compute(frame)

        act.post_process()
        return get_data_from_dtype(dslist, dtype=dtype)

    else:
        raise ValueError("Must provide command, mask, or indices")


def _distance_to_ref_or_point(traj=None,
                               mask=None,
                               ref=None,
                               point=None,
                               dtype='ndarray',
                               top=None,
                               frame_indices=None):
    """distance from atom (mask) to a given reference trajectory or a given point

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom selection
    ref : {None, tuple of two arrays, trajectory-like}
        if ref is not None and point is None: calculate distance from `mask` to
        each atom in ref.
        if ref is not None and point is not None: point will be used as
        vector direction
    point: array-like or tuple or list
        only use when ref is None
    dtype : return data type
    top : Topology
    frame_indices : array-like

    Returns
    -------
    out : ndarray, shape depends on input

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> pt.distance(traj, mask='@CA', ref=traj[0], dtype='ndarray').shape
    (10, 1, 12)
    >>> pt.distance(traj, mask='@CA', point=(0, 0, 0), dtype='ndarray').shape
    (10, 12)
    """

    top = top if top is not None else traj.top.copy()
    atom_indices = top.select(mask)

    if ref is not None:

        if hasattr(ref, 'xyz'):
            # single frame
            if ref.n_atoms == traj.n_atoms:
                # same system
                ref_xyz = ref.xyz[atom_indices]
                distances = np.empty((traj.n_frames, len(atom_indices)))

                for i, frame in enumerate(traj):
                    atom_xyz = frame.xyz[atom_indices]
                    diff = atom_xyz - ref_xyz
                    distances[i] = np.sqrt(np.sum(diff**2, axis=1))

                return distances
            else:
                # different systems - calculate distance to each atom in ref
                ref_xyz = ref.xyz
                distances = np.empty((traj.n_frames, len(atom_indices), ref.n_atoms))

                for i, frame in enumerate(traj):
                    atom_xyz = frame.xyz[atom_indices]
                    for j, atom_pos in enumerate(atom_xyz):
                        diff = ref_xyz - atom_pos
                        distances[i, j] = np.sqrt(np.sum(diff**2, axis=1))

                return distances
        else:
            raise ValueError("ref must be frame-like object with xyz attribute")

    elif point is not None:
        point = np.asarray(point, dtype=float)
        distances = np.empty((traj.n_frames, len(atom_indices)))

        for i, frame in enumerate(traj):
            atom_xyz = frame.xyz[atom_indices]
            diff = atom_xyz - point
            distances[i] = np.sqrt(np.sum(diff**2, axis=1))

        return distances
    else:
        raise ValueError("Either ref or point must be provided")


def pairwise_distance(traj=None,
                      mask='',
                      dtype='ndarray',
                      top=None,
                      frame_indices=None):
    """calculate pairwise distance from a number of atoms

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        mask to select atoms
    dtype : str, default 'ndarray'
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape (n_frames, n_atoms*(n_atoms-1)/2)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.pairwise_distance(traj(0, 2), mask=':3')
    >>> data.shape # distance of all atom pairs in residue 3
    (2, 66)
    """

    c_dslist = CpptrajDatasetList()
    c_action = c_action.Action_Pairwise()

    command = mask
    c_action.read_input(command, top=traj.top, dslist=c_dslist)
    c_action.setup(traj.top)

    for frame in traj:
        c_action.compute(frame)

    c_action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


def _check_command_type(command):
    """
    check command is string or array-like
    """
    if isinstance(command, str):
        return [command]
    else:
        return command


def _calculate_angles_for_int_array(traj, integer_array, n_frames, dtype):
    """
    Parameters
    ----------
    traj : trajectory-like
    integer_array : list or array, shape (n_triplets, 3)
        atom indices for angle calculation
    n_frames : int
    dtype : str

    Returns
    -------
    angles : array, shape (n_angles, n_frames)
    """
    angle_vals = np.empty((len(integer_array), n_frames))
    for i, (idx0, idx1, idx2) in enumerate(integer_array):
        angle_vals[i] = vector.angle(traj.xyz[:, idx0],
                                   traj.xyz[:, idx1],
                                   traj.xyz[:, idx2])
    return angle_vals


def _create_and_compute_action_list(list_of_commands: List[str],
                                  traj,
                                  action,
                                  dtype='ndarray'):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # need to add more examples
    '''
    from ..analysis.c_action.actionlist import ActionList

    actionlist = ActionList(list_of_commands, traj, action)
    with capture_stdout():
        actionlist.compute(traj)
    return get_data_from_dtype(actionlist.data, dtype=dtype)


def angle(traj=None,
          command='',
          mask=None,
          indices=None,
          dtype='ndarray',
          top=None,
          frame_indices=None):
    """compute angle

    Parameters
    ----------
    traj : Trajectory-like
    command : str or array-like
        cpptraj command or a list of cpptraj commands
    mask : str or array-like
        atom mask
    indices : array-like, shape (3,) or shape (n_angles, 3)
        atom index or a list of atom indices
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : {None, array-like}, optional

    Returns
    -------
    out : ndarray, shape depends on input

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # single command
    >>> data = pt.angle(traj, command="angle :3 :7 :8")
    >>> data.shape
    (10,)

    >>> # multiple commands
    >>> data = pt.angle(traj, command=["angle :3 :7 :8", "angle :4 :8 :9"])
    >>> data.shape
    (2, 10)

    >>> # using mask
    >>> data = pt.angle(traj, mask=":3 :7 :8")
    >>> data.shape
    (10,)

    >>> # using indices
    >>> data = pt.angle(traj, indices=[2, 6, 7])
    >>> data.shape
    (10,)
    """

    list_of_commands = _check_command_type(command)
    if len(list_of_commands) > 1 and not isinstance(command, str):
        # multiple commands
        return _create_and_compute_action_list(list_of_commands, traj,
                                             c_action.Action_Angle, dtype=dtype)

    if indices is not None:
        integer_array = np.asarray(indices, dtype=int)
        if integer_array.ndim == 1:
            if len(integer_array) != 3:
                raise ValueError("indices should have 3 items for angle")
            integer_array = integer_array.reshape(1, -1)

        angle_vals = _calculate_angles_for_int_array(traj, integer_array,
                                                   traj.n_frames, dtype)
        if len(integer_array) == 1:
            angle_vals = angle_vals.flatten()

        return angle_vals

    elif mask is not None:
        integer_array = get_list_of_commands(mask, top=traj.top)
        angle_vals = _calculate_angles_for_int_array(traj, integer_array,
                                                   traj.n_frames, dtype)
        if len(integer_array) == 1:
            angle_vals = angle_vals.flatten()

        return angle_vals

    else:
        # use command
        c_dslist = CpptrajDatasetList()
        action = c_action.Action_Angle()

        if command == '':
            raise ValueError("command can't be empty")

        action.read_input(command, top=traj.top, dslist=c_dslist)
        action.setup(traj.top)

        for frame in traj:
            action.compute(frame)

        action.post_process()
        return get_data_from_dtype(c_dslist, dtype=dtype)


def _dihedral_res(traj, mask=(), resid=0, dtype='ndarray', top=None):
    """used internally

    Return
    ------
    out : ndarray, shape (len(mask), traj.n_frames)

    Notes
    -----
    this method is used when `mask` has specific value
    User should not use this method. It will be changed or removed without
    announcement.
    """
    # use list instead of tuple so the order does not change
    dih_types = ['phi', 'psi', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5',
                'omega', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'nu1', 'nu2']
    arr = np.empty((len(mask), traj.n_frames))

    c_dslist = CpptrajDatasetList()

    for i, m in enumerate(mask):
        if m in dih_types:
            if resid != 0:
                command = f"dihedral :{resid}@{m}"
            else:
                command = f"dihedral {m}"
        else:
            command = f"dihedral {m}"

        c_action = c_action.Action_Dihedral()
        c_action.read_input(command, top=traj.top, dslist=c_dslist)
        c_action.setup(traj.top)

        for frame in traj:
            c_action.compute(frame)

        c_action.post_process()
        arr[i] = c_dslist[-1].values

    return arr


def _calculate_dihedrals_for_int_array(traj, integer_array, n_frames, dtype):
    """Calculate dihedrals for integer array indices

    Parameters
    ----------
    traj : trajectory-like
    integer_array : list or array, shape (n_dihedrals, 4)
        atom indices for dihedral calculation
    n_frames : int
    dtype : str

    Returns
    -------
    dihedrals : array, shape (n_dihedrals, n_frames)
    """
    dihedral_vals = np.empty((len(integer_array), n_frames))
    for i, (idx0, idx1, idx2, idx3) in enumerate(integer_array):
        dihedral_vals[i] = vector.dihedral(traj.xyz[:, idx0],
                                         traj.xyz[:, idx1],
                                         traj.xyz[:, idx2],
                                         traj.xyz[:, idx3])
    return dihedral_vals


def dihedral(traj=None,
             command='',
             mask=None,
             indices=None,
             resid=0,
             dtype='ndarray',
             top=None,
             frame_indices=None):
    """compute dihedral angle

    Parameters
    ----------
    traj : Trajectory-like
    command : str or List[str]
        cpptraj command
    mask : str, array-like, optional
        atom mask
    indices : array-like, shape (4,) or shape (n_dihedrals, 4)
        atom index
    resid : int, default 0
        residue id, start from 1. This is used to specify
        restraint for specific residue
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape depends on input

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # calc phi of 1st residue (index = 0)
    >>> data = pt.dihedral(traj, resid=1, mask='phi')
    >>> data.shape
    (10,)

    >>> # calc psi of 1st residue
    >>> data = pt.dihedral(traj, resid=1, mask='psi')
    >>> data.shape
    (10,)

    >>> # calc multiple dihedrals
    >>> data = pt.dihedral(traj, resid=1, mask=['phi', 'psi'])
    >>> data.shape
    (2, 10)

    >>> # use explicit command
    >>> data = pt.dihedral(traj, command='dihedral :1 :2 :3 :4')
    >>> data.shape
    (10,)
    """

    list_of_commands = _check_command_type(command)
    if len(list_of_commands) > 1 and not isinstance(command, str):
        # multiple commands
        return _create_and_compute_action_list(list_of_commands, traj,
                                             c_action.Action_Dihedral, dtype=dtype)

    if indices is not None:
        integer_array = np.asarray(indices, dtype=int)
        if integer_array.ndim == 1:
            if len(integer_array) != 4:
                raise ValueError("indices should have 4 items for dihedral")
            integer_array = integer_array.reshape(1, -1)

        dihedral_vals = _calculate_dihedrals_for_int_array(traj, integer_array,
                                                         traj.n_frames, dtype)
        if len(integer_array) == 1:
            dihedral_vals = dihedral_vals.flatten()

        return dihedral_vals

    elif mask is not None:
        if isinstance(mask, str):
            mask = [mask]

        if len(mask) == 1 and resid == 0:
            # single dihedral with mask
            integer_array = get_list_of_commands(mask[0], top=traj.top)
            dihedral_vals = _calculate_dihedrals_for_int_array(traj, integer_array,
                                                             traj.n_frames, dtype)
            return dihedral_vals.flatten()
        else:
            # use special handling for dihedral types
            return _dihedral_res(traj, mask, resid, dtype, top)
    else:
        # use command
        c_dslist = CpptrajDatasetList()
        action = c_action.Action_Dihedral()

        if command == '':
            raise ValueError("command can't be empty")

        action.read_input(command, top=traj.top, dslist=c_dslist)
        action.setup(traj.top)

        for frame in traj:
            action.compute(frame)

        action.post_process()
        return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def mindist(traj=None,
            mask='',
            dtype='ndarray',
            top=None,
            frame_indices=None):
    """compute mindist

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
    dtype : str, default 'ndarray'
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape (n_frames,)
    """
    c_action = c_action.Action_MinDist()
    c_dslist = CpptrajDatasetList()
    c_action.read_input(mask, top=traj.top, dslist=c_dslist)
    c_action.setup(traj.top)

    for frame in traj:
        c_action.compute(frame)

    c_action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def dihedral_rms(traj=None, mask="", dtype='ndarray', top=None, frame_indices=None, ref=None, extra_options=""):
    """Compute RMS of dihedral angles.

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, atom mask
    dtype : str, default 'ndarray'
        Output data type.
    top : Topology, optional
    frame_indices : array-like, optional
    ref : int or Frame, optional
        Reference frame for the calculation.
    extra_options : str, optional
        Additional cpptraj options for future extensibility.

    Returns
    -------
    DatasetList or ndarray
    """
    command_builder = CommandBuilder().add(mask)

    ref_name = None
    if ref is not None:
        ref_name = "myref"
        command_builder.add(f"ref {ref_name}")

    if extra_options:
        command_builder.add(extra_options)

    command = command_builder.build()

    action_datasets = CpptrajDatasetList()

    if ref is not None:
        ref_frame = get_reference(traj, ref)
        ref_dataset = action_datasets.add(DatasetType.REFERENCE_FRAME, name=ref_name)
        ref_dataset.top = ref_frame.top or traj.top
        ref_dataset.add_frame(ref_frame)

    action_datasets, _ = do_action(traj, command, c_action.Action_DihedralRMS, dslist=action_datasets)

    if ref is not None:
        action_datasets._pop(0)

    return get_data_from_dtype(action_datasets, dtype=dtype)