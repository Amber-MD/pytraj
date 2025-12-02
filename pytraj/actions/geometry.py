"""
Geometric analysis functions: distances, angles, dihedrals
"""
from .base import *

__all__ = [
    'distance', 'pairwise_distance', 'angle', 'dihedral', 'mindist',
    'dihedral_rms', 'distance_to_point', 'distance_to_reference', '_calculate_distance',
    '_calculate_angles_for_int_array', '_calculate_dihedrals_for_int_array', '_dihedral_res'
]


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
    """compute distance between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or a list of string or a 2D array-like of integers.
        If `mask` is a 2D-array, the `image` option is always `False`.
        In this case, make sure to `autoimage` your trajectory before
        calling `distance`.
    frame_indices : array-like, optional, default None
    dtype : return type, default 'ndarray'
    top : Topology, optional
    image : bool, default False
    n_frames : int, optional, default None
        only need to provide n_frames if ``traj`` does not have this info

    Returns
    -------
    1D ndarray if mask is a string
    2D ndarray, shape (n_atom_pairs, n_frames) if mask is a list of strings or an array

    Notes
    -----
    Be careful with Topology. If your topology has Box info but your traj does not, you
    would get weird output ([0.0, ...]). Make sure to use `image=False` in this method or
    set_nobox for Topology.


    Examples
    --------
    >>> import pytraj as pt
    >>> # calculate distance for two atoms, using amber mask
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> dist = pt.distance(traj, '@1 @3')

    >>> # calculate distance for two groups of atoms, using amber mask
    >>> dist = pt.distance(traj, '@1,37,8 @2,4,6')

    >>> # calculate distance between two residues, using amber mask
    >>> dist = pt.distance(traj, ':1 :10')

    >>> # calculate multiple distances between two residues, using amber mask
    >>> # distance between residue 1 and 10, distance between residue 3 and 20
    >>> # (when using atom string mask, index starts from 1)
    >>> dist = pt.distance(traj, [':1 :10', ':3 :20'])

    >>> # calculate distance for a series of atoms, using array for atom mask
    >>> # distance between atom 1 and 5, distance between atom 4 and 10 (index starts from 0)
    >>> dist = pt.distance(traj, [[1, 5], [4, 10]])
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
                      mask_1='',
                      mask_2='',
                      top=None,
                      dtype='ndarray',
                      frame_indices=None):
    '''compute pairwise distance between atoms in mask_1 and atoms in mask_2

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask_1: string or 1D array-like
    mask_2: string or 1D array-like
    ...

    Returns
    -------
    out_1 : numpy array, shape (n_frames, n_atom_1, n_atom_2)
    out_2 : atom pairs, shape=(n_atom_1, n_atom_2, 2)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> mat = pt.pairwise_distance(traj, '@CA', '@CB')

    Notes
    -----
    This method is only fast for small number of atoms.
    '''
    from itertools import product

    top_ = get_topology(traj, top)
    indices_1 = top_.select(mask_1) if isinstance(mask_1,
                                                  str) else mask_1
    indices_2 = top_.select(mask_2) if isinstance(mask_2,
                                                  str) else mask_2
    arr = np.array(list(product(indices_1, indices_2)))
    mat = distance(
        traj, mask=arr, dtype=dtype, top=top_, frame_indices=frame_indices)
    mat = mat.T
    return (mat.reshape(mat.shape[0], len(indices_1), len(indices_2)),
            arr.reshape(len(indices_1), len(indices_2), 2))


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





def dihedral(traj=None,
             mask="",
             top=None,
             dtype='ndarray',
             frame_indices=None,
             *args,
             **kwargs):
    """compute dihedral angle between two maskes
    ...
    """
    ensure_not_none_or_string(traj)
    command = mask

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_type = _check_command_type(command)

    if command_type == CommandType.STR:
        try:
            del kwargs['n_frames']
        except KeyError:
            pass
        action_datasets = CpptrajDatasetList()
        action = c_action.Action_Dihedral()
        action(command, traj, top=top, dslist=action_datasets, *args, **kwargs)
        return get_data_from_dtype(action_datasets, dtype)
    elif command_type == CommandType.LIST:
        return _create_and_compute_action_list(command, top, traj,
                                               c_action.Action_Dihedral, dtype,
                                               args, kwargs)
    elif command_type == CommandType.INT:
        integer_array = np.asarray(command)
        n_frames = kwargs.get('n_frames')
        n_frames = traj.n_frames if n_frames is None else n_frames
        return _calculate_dihedrals_for_int_array(traj, integer_array, n_frames, dtype)


@super_dispatch()
def mindist(traj=None,
            command="",
            top=None,
            dtype='ndarray',
            frame_indices=None):
    '''compute mindist

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.mindist(traj, '@CA @H')
    '''
    if not isinstance(command, str):
        command = array2d_to_cpptraj_maskgroup(command)
    command = "mindist " + command

    action_datasets, _ = do_action(traj, command, c_action.Action_NativeContacts)
    return get_data_from_dtype(action_datasets, dtype=dtype)[-1]


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