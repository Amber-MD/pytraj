"""
Topology manipulation functions: centering, alignment, imaging, etc.
"""
from .base import *
# Ensure _assert_mutable is available
from .base import _assert_mutable

__all__ = [
    'center_of_mass', 'center_of_geometry', 'align', 'align_principal_axis',
    'principal_axes', 'translate', 'rotate', 'autoimage', 'image', 'center',
    'strip', 'replicate_cell', 'rotate_dihedral', 'set_dihedral', 'set_velocity',
    'randomize_ions', 'fiximagedbonds', 'closest', 'atom_map', 'check_chirality',
    'scale', '_closest_iter'
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


@super_dispatch()
def center_of_mass(traj=None,
                   mask='*',
                   dtype='ndarray',
                   top=None,
                   frame_indices=None):
    """compute center of mass

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*' (all atoms)
        atom mask
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    output : ndarray, shape (n_frames, 3)
        center of mass coordinates
    """
    if dtype == 'ndarray':
        com_data = _calc_vector_center(traj, mask=mask, mass=True, top=top,
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


@super_dispatch()
def center_of_geometry(traj=None,
                       mask='*',
                       dtype='ndarray',
                       top=None,
                       frame_indices=None):
    """compute center of geometry (center of selected atoms coordinates)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*' (all atoms)
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
    if dtype == 'ndarray':
        com_data = _calc_vector_center(traj, mask=mask, mass=False, top=top,
                                       frame_indices=frame_indices)
        return com_data.reshape(traj.n_frames, 3)
    else:
        command = f"{mask} origin"
        dslist = CpptrajDatasetList()
        act = c_action.Action_Vector()
        act.read_input(command, top=traj.top, dslist=dslist)
        act.setup(traj.top)

        for frame in traj:
            act.compute(frame)

        act.post_process()
        return get_data_from_dtype(dslist, dtype=dtype)


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

        command = ' '.join((reference_command, mask_str, ref_mask_str, mass_str))

        if reference_topology is None:
            reference_topology = traj.top

        action_datasets = CpptrajDatasetList()
        action_datasets.add(DatasetType.REFERENCE, name=reference_name)
        action_datasets[0].top = reference_topology
        action_datasets[0].add_frame(ref)

        align_action = c_action.Action_Align()
        align_action.read_input(command, top=top, dslist=action_datasets)
        align_action.setup(top)

        for frame in traj:
            align_action.compute(frame)
        align_action.post_process()

        # remove ref
        action_datasets._pop(0)

        return traj


@super_dispatch()
def align_principal_axis(traj=None,
                         mask='',
                         dtype='ndarray',
                         top=None,
                         frame_indices=None):
    """align trajectory along principal axis

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    output : aligned trajectory
    """
    mut_traj = _assert_mutable(traj)

    for frame in mut_traj:
        pa = principal_axes(frame, mask, dorotation=True)
        frame.xyz = pa


def principal_axes(traj=None, mask='*', dorotation=False, mass=True, top=None):
    """compute principal axes

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*'
        atom mask
    dorotation : bool, default False
        if True, return array of aligned coordinates
        if False, return DatasetList storing principal axis vectors
    mass : bool, default True
        if True, use mass weighting
    top : Topology, optional

    Returns
    -------
    output : ndarray or DatasetList
    """

    if dorotation:
        # return aligned coordinate arrays
        mutraj = _assert_mutable(traj)
        command = mask

        if mass:
            command += " mass dorotation"
        else:
            command += " dorotation"

        action = c_action.Action_Principal()
        action.read_input(command, top=mutraj.top)
        action.setup(mutraj.top)

        for frame in mutraj:
            action.compute(frame)

        return mutraj.xyz
    else:
        # return dataset with principal axis vectors
        command = mask
        if mass:
            command += " mass"

        dslist = CpptrajDatasetList()
        action = c_action.Action_Principal()
        action.read_input(command, top=traj.top, dslist=dslist)
        action.setup(traj.top)

        for frame in traj:
            action.compute(frame)

        action.post_process()
        return get_data_from_dtype(dslist, dtype='dataset')


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
        c_dslist = CpptrajDatasetList()
        action = c_action.Action_Closest()
        action.read_input(command, top=traj.top, dslist=c_dslist)
        action.setup(traj.top)

        _closest_iter(c_action, traj)
        return get_data_from_dtype(c_dslist, dtype=dtype)


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
    mut_traj = _assert_mutable(traj)

    action = c_action.Action_Translate()
    action.read_input(command, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


@super_dispatch()
def scale(traj=None, command="", frame_indices=None, top=None):
    """do coordinate scaling"""
    mut_traj = _assert_mutable(traj)

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
    mut_traj = _assert_mutable(traj)

    action = c_action.Action_Rotate()
    action.read_input(command, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


@super_dispatch()
def autoimage(traj, mask="", frame_indices=None, top=None):
    """autoimage for trajectories

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
    >>> _ = pt.autoimage(traj)
    """
    mut_traj = _assert_mutable(traj)

    action = c_action.Action_AutoImage()
    action.read_input(mask, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


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
    mut_traj = _assert_mutable(traj)

    action = c_action.Action_Image()
    action.read_input(mask, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


@super_dispatch()
def center(traj=None,
           mask='origin',
           mass=False,
           frame_indices=None,
           top=None):
    """center atoms

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default 'origin'
        cpptraj command
    mass : bool, default False
        if True, use mass-weighted centering
    frame_indices : array-like, optional
    top : Topology, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> xyz_old = traj.xyz[0].copy()
    >>> _ = pt.center(traj, mask='origin')
    >>> (traj.xyz[0] == xyz_old).all()
    False
    """
    mut_traj = _assert_mutable(traj)

    command = mask
    if mass:
        command += ' mass'

    action = c_action.Action_Center()
    action.read_input(command, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


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
    if hasattr(obj, 'top'):
        # Trajectory-like
        indices = obj.top.select(mask, invert=True)
        return obj[indices]
    else:
        # Topology
        indices = obj.select(mask, invert=True)
        return obj[indices]


@super_dispatch()
def randomize_ions(traj,
                   mask='',
                   around_mask='',
                   min_dist=5.0,
                   options='',
                   dtype='ndarray',
                   frame_indices=None,
                   top=None):
    """randomize ions coordinates

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    around_mask : str, optional
        around mask
    min_dist : float, default 5.0
        minimum distance
    options : str, optional
        extra options
    dtype : str, default 'ndarray'
        return data type
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    out : Trajectory
    """
    command = mask
    if around_mask:
        command += f" around {around_mask}"
    if min_dist != 5.0:
        command += f" by {min_dist}"
    if options:
        command += f" {options}"

    mut_traj = _assert_mutable(traj)

    action = c_action.Action_RandomizeIons()
    action.read_input(command, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


@super_dispatch()
def replicate_cell(traj=None,
                   replicates=(1, 1, 1),
                   mask='*',
                   dtype='ndarray',
                   frame_indices=None,
                   top=None):
    """replicate unit cell

    Parameters
    ----------
    traj : Trajectory-like
    replicates : tuple of ints, default (1, 1, 1)
        how many times to replicate along each dimension
    mask : str, default '*'
        atom mask
    dtype : str, default 'ndarray'
        return data type
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    out : Trajectory
    """
    nx, ny, nz = replicates
    command = f"{mask} {nx}x{ny}x{nz}"

    mut_traj = _assert_mutable(traj)

    action = c_action.Action_ReplicateCell()
    action.read_input(command, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    action.post_process()

    return mut_traj


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
    """set specific dihedral angle to specific value"""

    if dihedral_type is None:
        dihedral_type = 'phi'

    mut_traj = _assert_mutable(traj)

    if resid != 0:
        mask = ":{:d}@{:s}".format(resid, dihedral_type)
    else:
        mask = dihedral_type

    mask += " {}".format(deg)

    action = c_action.Action_SetDihedral()
    action.read_input(mask, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


def set_velocity(traj, temperature=298, ig=10, options=''):
    """assign Maxwell-Boltzmann velocities

    Parameters
    ----------
    traj : Trajectory-like
    temperature : float, default 298
        temperature in K
    ig : int, default 10
        random seed for Maxwell-Boltzmann distribution
    options : str, optional
        extra options

    Returns
    -------
    traj : Trajectory with velocities assigned

    Notes
    -----
    Need velocity information
    """
    mut_traj = _assert_mutable(traj)

    command = f"ig {ig} temp0 {temperature} " + options

    action = c_action.Action_SetVelocity()
    action.read_input(command, top=mut_traj.top)
    action.setup(mut_traj.top)

    for frame in mut_traj:
        action.compute(frame)

    return mut_traj


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
    mut_traj = _assert_mutable(traj)

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

    target = dataset_list.add('reference', name='my_target')
    target.top = traj.top
    target.append(traj[0])

    refset = dataset_list.add('reference', name='my_ref')
    refset.top = ref.top if ref.top is not None else traj.top
    if not isinstance(ref, Frame):
        ref_frame = ref[0]
    else:
        ref_frame = ref
    refset.append(ref_frame)

    with capture_stdout() as (out, err):
        act(command, traj, top=traj.top, dslist=dataset_list)
    act.post_process()

    # free memory of two reference
    dataset_list._pop(0)
    dataset_list._pop(0)

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
        return Trajectory(
            xyz=np.array([frame.xyz.copy() for frame in fiter]),
            top=new_top.copy())
    else:
        # iterator
        return (fiter, new_top.copy())