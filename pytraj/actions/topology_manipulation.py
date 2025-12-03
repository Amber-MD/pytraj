"""
Topology manipulation functions: centering, alignment, imaging, etc.
"""
from .base import *
from .base import add_reference_dataset
# Ensure _assert_mutable is available
from .base import _assert_mutable, _ensure_mutable
from ..builder.build import make_structure

__all__ = [
    'center_of_mass', 'center_of_geometry', 'align', 'align_principal_axis',
    'principal_axes', 'translate', 'rotate', 'autoimage', 'image', 'center',
    'strip', 'replicate_cell', 'rotate_dihedral', 'set_dihedral',
    'set_velocity', 'randomize_ions', 'fiximagedbonds', 'closest', 'atom_map',
    'check_chirality', 'scale', '_closest_iter'
]


def _calc_vector_center(traj=None,
                        mask='*',
                        mass=True,
                        top=None,
                        frame_indices=None):
    """compute 'center' of selected atoms."""
    command = "center " + mask
    if mass:
        command += " mass"

    dslist = CpptrajDatasetList()
    act = c_action.Action_Vector()
    # Ensure top is a valid Topology object
    if top is None:
        top = traj.top
    act.read_input(command, top=top, dslist=dslist)
    act.setup(top)

    for frame in traj:
        act.compute(frame)

    act.post_process()
    return dslist[0].values


def center_of_mass(traj=None,
                   mask='',
                   top=None,
                   dtype='ndarray',
                   frame_indices=None):
    """compute center of mass

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '' (all atoms)
        atom mask
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    output : ndarray, shape (n_frames, 3)
        center of mass coordinates

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # compute center of mass residue 3 for first 2 frames.
    array([[-0.661702  ,  6.69124347,  3.35159413],
           [ 0.5620708 ,  7.82263042, -0.72707798]])
    """
    # note: do not use super_dispatch for this method since
    # we already use for _calc_vector_center
    if dtype == 'ndarray':
        com_data = _calc_vector_center(traj=traj,
                                       mask=mask,
                                       top=top,
                                       mass=True,
                                       frame_indices=frame_indices)
        return com_data.reshape(traj.n_frames, 3)
    else:
        command = "center " + mask + " mass"
        dslist = CpptrajDatasetList()
        act = c_action.Action_Vector()
        if top is None:
            top = traj.top
        act.read_input(command, top=top, dslist=dslist)
        act.setup(top)

        for frame in traj:
            act.compute(frame)

        act.post_process()
        return get_data_from_dtype(dslist, dtype=dtype)


@register_pmap
@super_dispatch()
def center_of_geometry(traj=None,
                       mask="",
                       top=None,
                       dtype='ndarray',
                       frame_indices=None):
    """compute center of geometry (center of selected atoms coordinates)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '' (all atoms)
        atom mask
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    output : ndarray, shape (n_frames, 3)
        center of geometry coordinates
    """
    atom_mask_obj = top(mask)
    action_datasets = CpptrajDatasetList()
    action_datasets.add(DatasetType.VECTOR)

    for frame in iterframe_master(traj):
        action_datasets[0].append(frame.center_of_geometry(atom_mask_obj))
    return get_data_from_dtype(action_datasets, dtype=dtype)


@super_dispatch()
def align(traj,
          mask='',
          ref=0,
          ref_mask='',
          mass=False,
          top=None,
          frame_indices=None):
    """align (superpose) trajectory to given reference

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '' (all atoms)
    ref : {int, Frame}, default 0 (first frame)
    ref_mask : str, default ''
        if not given, use traj's mask
        if given, use it
    mass : Bool, default False
        if True, mass-weighted
        if False, no mas-weighted
    frame_indices : {None, array-like}, default None
       if given, only compute RMSD for those

    Examples
    --------

    Notes
    -----
    versionadded: 1.0.6
    """
    if isinstance(traj, TrajectoryIterator):
        return traj.superpose(mask=mask, ref=ref, ref_mask=ref_mask, mass=mass)
    else:
        mask_str = mask
        ref_mask_str = ref_mask
        reference_topology = ref.top if hasattr(ref, 'top') else top
        mass_str = 'mass' if mass else ''

        reference_name = 'myref'
        reference_command = 'ref {}'.format(reference_name)

        command = ' '.join(
            (reference_command, mask_str, ref_mask_str, mass_str))

        if reference_topology is None:
            reference_topology = traj.top

        action_datasets = CpptrajDatasetList()
        add_reference_dataset(action_datasets, reference_name, ref, reference_topology)

        align_action = c_action.Action_Align()
        align_action.read_input(command, top=top, dslist=action_datasets)
        align_action.setup(top)

        for frame in traj:
            align_action.compute(frame)
        align_action.post_process()

        # remove ref
        action_datasets.remove_at(0)

        return traj


def align_principal_axis(traj=None,
                         mask="*",
                         top=None,
                         frame_indices=None,
                         mass=False):
    # TODO : does not match with cpptraj output
    # rmsd_nofit ~ 0.5 for md1_prod.Tc5b.x, 1st frame
    """
    Notes
    -----
    apply for mutatble traj (Trajectory, Frame)
    """
    mut_traj = _ensure_mutable(traj)
    mass_ = 'mass' if mass else ''
    command = ' '.join((mask, " dorotation", mass_))

    # Handle frame_indices manually since do_action doesn't support it
    c_dslist = CpptrajDatasetList()
    topology = get_topology(mut_traj, top)
    action = c_action.Action_Principal()
    action.read_input(command, top=topology)
    action.setup(topology)

    # Process only specified frames or all frames
    if frame_indices is not None:
        for idx, frame in enumerate(mut_traj):
            if idx in frame_indices:
                action.compute(frame)
    else:
        for frame in mut_traj:
            action.compute(frame)

    action.post_process()
    return mut_traj


def principal_axes(traj=None, mask='*', dorotation=False, mass=True, top=None):
    # TODO: update doc please
    """compute principal axes

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*' (all atoms)
    mass: bool, defaul True
    if `dorotation`, the system will be aligned along principal axes (apply for mutable system)
    top : Topology, optional

    Returns
    -------
    out_0 : numpy array, shape=(n_frames, 3, 3)
    out_1: numpy array with shape=(n_frames, 3)
    """
    command_elements = [mask]
    if 'name' not in mask:
        command_elements.append('name pout')
    if dorotation:
        command_elements.append('dorotation')
    if mass:
        command_elements.append('mass')

    command = ' '.join(command_elements)
    action_datasets, _ = do_action(traj, command, c_action.Action_Principal)

    principal_axes = action_datasets[0].values
    principal_values = action_datasets[1].values

    return principal_axes, principal_values


def _closest_iter(act, traj):
    """internal use"""
    for frame in traj:
        act.compute(frame)

    act.post_process()


@super_dispatch()
def closest(traj=None,
            mask='',
            top=None,
            solvent_mask=':WAT,Na+,Cl-,K+',
            n_solvents=10,
            closestout='',
            dtype='trajectory',
            frame_indices=None):
    """find closest water to solute

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        solute mask, default ''
    top : Topology, optional
    solvent_mask : str, default ':WAT,Na+,Cl-,K+'
        solvent mask
    n_solvents : int, default 10
        number of closest solvents
    closestout : str, optional
        output filename
    dtype : str, default 'trajectory'
        return data type
    frame_indices : array-like, optional

    Returns
    -------
    traj_out : Trajectory
        if dtype is 'trajectory', return stripped trajectory
        with closest waters. Otherwise return DatasetList
    """
    command = f"{n_solvents} {solvent_mask}"
    if mask:
        command = f"{mask} {command}"
    if closestout:
        command += f" closestout {closestout}"

    if dtype == 'trajectory':
        traj_mut = traj[:]
        action = c_action.Action_Closest()
        action.read_input(command, top=traj_mut.top)
        action.setup(traj_mut.top)

        _closest_iter(c_action, traj_mut)

        return traj_mut
    else:
        action_datasets, _ = do_action(traj, command, c_action.Action_Closest)
        return get_data_from_dtype(action_datasets, dtype=dtype)


@super_dispatch()
def translate(traj=None, command="", top=None, frame_indices=None):
    """translate atoms by given vector

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        cpptraj command, example: 'x 100, y 100, z 0.0'
    top : Topology, optional
    frame_indices : array-like, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> xyz_old = traj.xyz[0, 0].copy()
    >>> _ = pt.translate(traj, "x 1.2")
    >>> (traj.xyz[0, 0] - xyz_old) # doctest: +SKIP
    array([ 1.20000005,  0.        ,  0.        ], dtype=float32)
    """
    mut_traj = _ensure_mutable(traj)

    action = c_action.Action_Translate()
    action.read_input(command, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


@super_dispatch()
def scale(traj=None, command="", frame_indices=None, top=None):
    """do coordinate scaling"""
    mut_traj = _ensure_mutable(traj)

    action = c_action.Action_Scale()
    action.read_input(command, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


@super_dispatch()
def rotate(traj=None, command="", frame_indices=None, top=None):
    """rotate atoms

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        cpptraj rotate command
    frame_indices : array-like, optional
    top : Topology, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # rotate along x, y, z by 30 degree
    >>> _ = pt.rotate(traj, "x 30 y 30 z 30")
    """
    mut_traj = _ensure_mutable(traj)

    action = c_action.Action_Rotate()
    action.read_input(command, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


def autoimage(traj, mask="", frame_indices=None, top=None):
    '''perform autoimage and return the updated-coordinate traj

    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()[:]
    >>> traj = pt.autoimage(traj)
    '''
    top = top or traj.top
    assert top.has_box(), "Topology must have box information"
    _assert_mutable(traj)
    command = mask
    do_action(traj, command, c_action.Action_AutoImage, top=top)
    return traj


@super_dispatch()
def image(traj, mask="", frame_indices=None, top=None):
    """move atoms to primary unit cell

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    frame_indices : array-like, optional
    top : Topology, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> _ = pt.image(traj, "origin")
    """
    mut_traj = _ensure_mutable(traj)

    action = c_action.Action_Image()
    action.read_input(mask, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


@super_dispatch()
def center(traj=None,
           mask="",
           center='box',
           mass=False,
           top=None,
           frame_indices=None):
    """Center coordinates in `mask` to specified point.

    Parameters
    ----------
    traj : Trajectory-like or Frame iterator
    mask : str, mask
    center : str, {'box', 'origin', array-like}, default 'box'
        if 'origin', center on coordinate origin (0, 0, 0)
        if 'box', center on box center
        if array-like, center on that point
    mass : bool, default: False
        if True, use mass weighted
    top : Topology, optional, default: None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # load all frames to memory so we can 'mutate' them
    >>> traj = traj[:]
    >>> # all atoms, center to box center (x/2, y/2, z/2)
    >>> traj = pt.center(traj)

    >>> # center at origin, use @CA
    >>> traj = pt.center(traj, '@CA', center='origin')

    >>> # center to box center, use mass weighted
    >>> traj = pt.center(traj, mass=True)
    >>> traj = pt.center(traj, ':1', mass=True)

    Returns
    -------
    updated traj

    See also
    --------
    pytraj.translate
    """
    valid_centers = ['box', 'origin']

    if isinstance(center, (list, tuple)):
        center = 'point ' + ' '.join(map(str, center))
    elif center.lower() not in valid_centers:
        raise ValueError(f'center must be one of {valid_centers}')

    center_option = '' if center == 'box' else center
    mass_option = 'mass' if mass else ''
    command = ' '.join((mask, center_option, mass_option))

    if isinstance(traj, TrajectoryIterator):
        return traj.center(command)

    action_center = c_action.Action_Center()
    action_center(command, traj, top=top)

    return traj


def strip(obj, mask):
    """return a new trajectory with atoms stripped

    Parameters
    ----------
    obj : Trajectory-like or Topology
    mask : str
        atom mask

    Returns
    -------
    out : Trajectory or Topology

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> traj.n_atoms
    5293
    >>> traj2 = pt.strip(traj, '@H=') # strip all Hydrogens
    >>> traj2.n_atoms
    2668
    >>> traj.n_atoms
    5293
    """
    if isinstance(obj, str) and not isinstance(mask, str):
        obj, mask = mask, obj

    kept_mask = '!(' + mask + ')'

    if isinstance(obj, (Topology, Trajectory)):
        # return new Topology or new Trajectory
        return obj[kept_mask]
    elif isinstance(obj, TrajectoryIterator):
        # return a FrameIterator
        return obj(mask=kept_mask)
    elif hasattr(obj, 'mask'):
        obj.mask = kept_mask
        return obj
    else:
        raise ValueError('object must be either Trajectory or Topology')


@super_dispatch()
def randomize_ions(traj,
                   mask,
                   around,
                   by,
                   overlap,
                   seed=1,
                   top=None,
                   frame_indices=None):
    """randomize_ions for given Frame with Topology

    Parameters
    ----------
    traj : Trajectory-like or a Frame
        ``traj`` must be mutable
    mask : str
        cpptraj command
    frame_indices : {None, array-like}, optional
    top : Topology, optional (only needed if ``traj`` does not have it)

    Examples
    --------
    >>> pt.randomize_ions(traj, mask='@Na+', around=':1-16', by=5.0, overlap=3.0, seed=113698) # doctest: +SKIP
    """
    _assert_mutable(traj)
    around_ = 'around ' + str(around)
    by_ = 'by ' + str(by)
    overlap_ = 'overlap ' + str(overlap)
    seed_ = 'seed ' + str(seed)
    command = ' '.join((mask, around_, by_, overlap_, seed_))

    do_action(traj, command, c_action.Action_RandomizeIons, top=top)
    return traj


@super_dispatch()
def replicate_cell(traj=None,
                   mask="",
                   direction='all',
                   frame_indices=None,
                   top=None):
    '''create a trajectory where the unit cell is replicated in 1 or more direction (up to 27)

    Parameters
    ----------
    traj : Trajectory-like or Frame iterator
    mask : str, default: ""
    direction : str or array-like, default: "all"
        how to replicate.

    Returns
    -------
    traj : pytraj.Trajectory

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> new_traj = pt.replicate_cell(traj, direction='all')
    >>> new_traj = pt.replicate_cell(traj, direction='dir 001 dir 111')
    >>> new_traj = pt.replicate_cell(traj, direction='dir 001 dir 1-10')
    >>> new_traj = pt.replicate_cell(traj, direction='dir 001 dir 1-10')

    >>> # similiar usage
    >>> new_traj = pt.replicate_cell(traj, direction=('001', '0-10'))
    >>> new_traj = pt.replicate_cell(traj, direction=['001', '0-10'])
    '''
    if isinstance(direction, str):
        formatted_direction = direction
    elif isinstance(direction, (list, tuple)):
        formatted_direction = 'dir ' + ' dir '.join(direction)
    else:
        raise ValueError('direction must be a string or list/tuple of strings')

    command = f'name tmp_cell {formatted_direction} {mask}'
    action_datasets, _ = do_action(traj, command, c_action.Action_ReplicateCell)

    return Trajectory(xyz=action_datasets[0].xyz, top=action_datasets[0].top)


def rotate_dihedral(traj=None, mask="", top=None):
    # change to pt.rotate_dihedral(traj, res=0,
    #              mask=("O4'", "C1'", "N9", "C4"), deg=120)?
    """
    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_rna()[:]
    >>> traj = pt.rotate_dihedral(traj, "3:chin:120") # rotate chin of res 3 to 120 deg
    >>> traj = pt.rotate_dihedral(traj, "1:O4':C1':N9:C4:120") # rotate dihedral with given mask

    Returns
    -------
    updated traj

    Notes
    -----
    Syntax and method's name might be changed
    """
    from ..utils.get_common_objects import get_topology

    _assert_mutable(traj)
    top_ = get_topology(traj, top)

    if "custom:" in mask:
        command = mask
    else:
        command = "custom:" + mask

    act = c_action.Action_MakeStructure()

    act(command, traj, top=top_)
    return traj


def set_dihedral(traj, resid=0, dihedral_type=None, deg=0, top=None):
    '''

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # make mutable traj by loading all frames to disk
    >>> traj = traj[:]
    >>> traj = pt.set_dihedral(traj, resid=2, dihedral_type='phi', deg=60)

    Returns
    -------
    updated traj
    '''
    if not isinstance(resid, str):
        resid = str(resid + 1)
    deg = str(deg)

    command = ':'.join((dihedral_type, resid, dihedral_type, deg))
    make_structure(traj, command)
    return traj


def set_velocity(traj, temperature=298, ig=10, options=''):
    """

    Notes
    -----
    Added in v2.0.1
    """
    command = "tempi {} ig {} {}".format(temperature, ig, options)

    top = traj.top

    act = c_action.Action_SetVelocity()
    act.read_input(command, top=top)
    act.setup(top)

    if traj.velocities is None:
        traj.velocities = np.empty(traj.xyz.shape)
    for index, frame in enumerate(traj):
        new_frame = act.compute(frame, get_new_frame=True)
        traj.velocities[index] = new_frame.velocity
    return traj


def fiximagedbonds(traj, mask=''):
    """fix imaged bonds

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask

    Returns
    -------
    traj : modified trajectory
    """
    mut_traj = _ensure_mutable(traj)

    action = c_action.Action_FixImagedBonds()
    action.read_input(mask, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


def atom_map(traj, ref, rmsfit=False):
    ''' Limited support for cpptraj atommap

    Parameters
    ----------
    traj : Trajectory-like
    ref : Trajectory-like with one frame
    rmsfit : bool, default False
        if True, compute rmsfit

    Notes
    -----
    This method in pytraj is not mature yet.

    Returns
    -------
    out : Tuple[str, np.ndarray]
        (mask_out, rmsd data if rmsfit=True)
    '''
    act = c_action.Action_AtomMap()
    options = 'rmsfit rmsout rmsout.dat' if rmsfit else ''
    command = ' '.join(('my_target my_ref', options))
    dataset_list = CpptrajDatasetList()

    add_reference_dataset(dataset_list, 'my_target', traj[0], traj.top)

    if not isinstance(ref, Frame):
        ref_frame = ref[0]
    else:
        ref_frame = ref
    add_reference_dataset(dataset_list, 'my_ref', ref_frame, ref.top if ref.top is not None else traj.top)

    with capture_stdout() as (out, err):
        act(command, traj, top=traj.top, dslist=dataset_list)
    act.post_process()

    # free memory of two reference
    dataset_list.remove_at(0)
    dataset_list.remove_at(0)

    return (out.read(), get_data_from_dtype(dataset_list, dtype='ndarray'))


def check_chirality(traj, mask='', dtype='dict'):
    """check chirality

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    dtype : str, default 'dict'
        return data type

    Returns
    -------
    out : dict or DatasetList
    """
    c_dslist = CpptrajDatasetList()
    action = c_action.Action_CheckChirality()
    action.read_input(mask, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


def _closest_iter(act, traj):
    '''Helper function for closest action

    Parameters
    ----------
    act : Action object
    traj : Trajectory-like
    '''
    for frame in iterframe_master(traj):
        new_frame = act.compute(frame, get_new_frame=True)
        yield new_frame


@register_openmp
@super_dispatch()
def closest(traj=None,
            mask='*',
            solvent_mask=None,
            n_solvents=10,
            frame_indices=None,
            dtype='iterator',
            top=None):
    """return either a new Trajectory or a frame iterator. Keep only ``n_solvents`` closest to mask

    Parameters
    ----------
    traj : Trajectory-like | list of Trajectory-like/frames | frame iterator | chunk iterator
    mask: str, default '*' (all solute atoms)
    top : Topology-like object, default=None, optional
    dtype : {'iterator', 'trajectory'}, default 'iterator'
        if 'iterator', return a tuple of Frame iterator and new Toplogy. 'iterator' is good for streaming large trajectory.
        if 'trajectory', return a new Trajectory.

    Returns
    -------
    out : (check above)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # obtain new traj, keeping only closest 100 waters
    >>> # to residues 1 to 13 (index starts from 1) by distance to the first atom of water
    >>> t = pt.closest(traj, mask='@CA', n_solvents=10)
    """
    # check if top has solvent
    c_dslist = CpptrajDatasetList()

    command = str(n_solvents) + ' ' + mask

    act = c_action.Action_Closest()

    if solvent_mask is not None:
        top = top.copy()
        top.set_solvent(solvent_mask)

    has_solvent = False
    for mol in top.mols:
        if mol.is_solvent():
            has_solvent = True
            break
    if not has_solvent:
        raise RuntimeError("Topology does not have solvent")

    act.read_input(command, top, dslist=c_dslist)
    new_top = act.setup(top, get_new_top=True)

    fiter = _closest_iter(act, traj)

    if dtype == 'trajectory':
        return Trajectory(xyz=np.array([frame.xyz.copy() for frame in fiter]),
                          top=new_top.copy())
    else:
        # iterator
        return (fiter, new_top.copy())
