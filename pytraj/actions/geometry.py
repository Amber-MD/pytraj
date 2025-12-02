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





def distance(traj=None,
             mask="",
             frame_indices=None,
             dtype='ndarray',
             top=None,
             image=False,
             n_frames=None):
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
    ensure_not_none_or_string(traj)
    command = mask

    traj = get_fiterator(traj, frame_indices)
    topology = get_topology(traj, top)
    noimage_str = 'noimage' if not image or (hasattr(traj, 'crdinfo') and not traj.crdinfo['has_box'] and topology.has_box()) else ''

    command_array = np.asarray(command)

    if isinstance(command, (list, tuple, str, np.ndarray, int)):
        if 'int' in command_array.dtype.name:
            integer_array = command_array
            frame_count = traj.n_frames if n_frames is None else n_frames
            return _calculate_distance(traj, integer_array, frame_count, dtype)
        else:
            command_list = get_list_of_commands(command) if not isinstance(command, np.ndarray) else command
            cpptraj_action_datasets = CpptrajDatasetList()
            action_list = ActionList()

            for command in command_list:
                if noimage_str:
                    command = ' '.join((command, noimage_str))
                action_list.add(c_action.Action_Distance(), command, topology, dslist=cpptraj_action_datasets)

            action_list.compute(traj)
            return get_data_from_dtype(cpptraj_action_datasets, dtype)
    else:
        raise ValueError(
            "command must be a string, a list/tuple of strings, or "
            "a numpy 2D array")


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
    action = c_action.Action_Pairwise()

    command = mask
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    c_action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


class CommandType(StrEnum):
    INT = 'int'
    STR = 'str'
    LIST = 'list'

def _check_command_type(command):
    command_array = np.asarray(command)
    if 'int' in command_array.dtype.name:
        return CommandType.INT
    elif isinstance(command, str):
        return CommandType.STR
    elif isinstance(command, (list, tuple, np.ndarray)):
        return CommandType.LIST
    else:
        raise ValueError("command must be a string, a list/tuple of strings, or a numpy 2D array")





def _create_and_compute_action_list(list_of_commands: List[str],
                                    top,
                                    traj,
                                    action: Callable,
                                    dtype: str,
                                    args: tuple,
                                    kwargs: dict):
    action_datasets = CpptrajDatasetList()
    action_list = ActionList()

    for command in list_of_commands:
        action_list.add(
            action(),
            command,
            top,
            dslist=action_datasets,
            *args,
            **kwargs)
    action_list.compute(traj)
    return get_data_from_dtype(action_datasets, dtype)
def angle(traj=None,
          mask="",
          frame_indices=None,
          dtype='ndarray',
          top=None,
          *args,
          **kwargs):
    """compute angle between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or array
    top : Topology, optional
    dtype : return type, defaul 'ndarray'

    Returns
    -------
    1D ndarray if mask is a string
    2D ndarray, shape (n_atom_pairs, n_frames) if mask is a list of strings or an array

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # calculate angle for three atoms, using amber mask
    >>> pt.angle(traj, '@1 @3 @10')
    array([  98.06193365,   95.75979717,  105.59626997,  107.64079091,
             94.93516228,  102.06028369,  109.3479057 ,  110.11532478,
            101.86104127,  110.48992512])

    >>> # calculate angle for three groups of atoms, using amber mask
    >>> angles = pt.angle(traj, '@1,37,8 @2,4,6 @5,20')

    >>> # calculate angle between three residues, using amber mask
    >>> angles = pt.angle(traj, ':1 :10 :20')

    >>> # calculate multiple angles between three residues, using amber mask
    >>> # angle between residue 1, 10, 20, angle between residue 3, 20, 30
    >>> # (when using atom string mask, index starts from 1)
    >>> angles = pt.angle(traj, [':1 :10 :20', ':3 :20 :30'])

    >>> # calculate angle for a series of atoms, using array for atom mask
    >>> # angle between atom 1, 5, 8, distance between atom 4, 10, 20 (index starts from 0)
    >>> angles = pt.angle(traj, [[1, 5, 8], [4, 10, 20]])
    """
    ensure_not_none_or_string(traj)
    command = mask

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_type = _check_command_type(command)

    if command_type == CommandType.STR:
        # need to remove 'n_frames' keyword since Action._master does not use
        # it
        try:
            del kwargs['n_frames']
        except KeyError:
            pass
        # cpptraj mask for action
        action_datasets = CpptrajDatasetList()
        action = c_action.Action_Angle()
        action(command, traj, top=top, dslist=action_datasets, *args, **kwargs)
        return get_data_from_dtype(action_datasets, dtype)
    elif command_type == CommandType.LIST:
        return _create_and_compute_action_list(command, top, traj,
                                            c_action.Action_Angle, dtype, args, kwargs)
    elif command_type == CommandType.INT:
        integer_array = np.asarray(command)
        n_frames = kwargs.get('n_frames')
        n_frames = traj.n_frames if n_frames is None else n_frames
        return _calculate_angles_for_int_array(traj, integer_array, n_frames, dtype)


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

        action = c_action.Action_Dihedral()
        action.read_input(command, top=traj.top, dslist=c_dslist)
        action.setup(traj.top)

        for frame in traj:
            action.compute(frame)

        action.post_process()
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
             range360=False,
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

        if range360:
            command += " range360"
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
    action = c_action.Action_MinDist()
    c_dslist = CpptrajDatasetList()
    action.read_input(mask, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

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