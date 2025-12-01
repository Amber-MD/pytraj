"""
Topology manipulation functions: centering, alignment, imaging, etc.
"""
from .base import *

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
    command = f"{mask} origin"
    if mass:
        command += " mass"

    dslist = CpptrajDatasetList()
    act = c_action.Action_Vector()
    act.read_input(command, top=traj.top, dslist=dslist)
    act.setup(traj.top)

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
        command = f"{mask} origin mass"
        dslist = CpptrajDatasetList()
        act = c_action.Action_Vector()
        act.read_input(command, top=traj.top, dslist=dslist)
        act.setup(traj.top)

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
          mask=None,
          ref=0,
          mass=False,
          rms_name='',
          frame_indices=None,
          top=None):
    """align trajectory to a reference frame

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        mask to select atoms for fitting
    ref : {int, Frame, Trajectory}, default 0 (first frame)
        Reference frame
    mass : bool, default False
        if True, use mass-weighted alignment
    rms_name : str, optional
        root mean squared name to be written to DatasetList
    frame_indices : array-like, optional
        frame indices
    top : Topology, optional

    Returns
    -------
    dslist : DatasetList

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> traj.xyz[0, 0]
    array([-1.88900006,  9.1590004 ,  0.24699999], dtype=float32)
    >>> _ = pt.align(traj, mask='!@H=', ref=traj[2])
    >>> traj.xyz[0, 0]  # coordinates changed
    array([ 6.97324038, -8.83906818,  5.22135925], dtype=float32)
    """
    from .rmsd import rmsd

    if mask is None:
        mask = "*"

    # make a copy of reference
    ref_frame = get_reference(traj, ref)

    command = 'rms2d'

    if rms_name:
        command = command + ' ' + rms_name
    if mass:
        command = command + ' mass'

    # create mutable trajectory
    mut_traj = _assert_mutable(traj)

    superpose(mut_traj, mask=mask, ref=ref_frame, mass=mass, frame_indices=frame_indices)

    dslist = rmsd(mut_traj, mask=mask, ref=ref_frame, mass=mass, dtype='dataset',
                  frame_indices=frame_indices)

    dslist.add(DatasetType.COORDS, name='_DEFAULTCRD_')
    dslist[-1].top = mut_traj.top
    for frame in mut_traj:
        dslist[-1].append(frame)

    return dslist


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
        c_action = c_action.Action_Closest()
        c_action.read_input(command, top=traj_mut.top)
        c_action.setup(traj_mut.top)

        _closest_iter(c_action, traj_mut)

        return traj_mut
    else:
        c_dslist = CpptrajDatasetList()
        c_action = c_action.Action_Closest()
        c_action.read_input(command, top=traj.top, dslist=c_dslist)
        c_action.setup(traj.top)

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

    c_action = c_action.Action_Translate()
    c_action.read_input(command, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

    return mut_traj


@super_dispatch()
def scale(traj=None, command="", frame_indices=None, top=None):
    """do coordinate scaling"""
    mut_traj = _assert_mutable(traj)

    c_action = c_action.Action_Scale()
    c_action.read_input(command, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

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

    c_action = c_action.Action_Rotate()
    c_action.read_input(command, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

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

    c_action = c_action.Action_AutoImage()
    c_action.read_input(mask, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

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

    c_action = c_action.Action_Image()
    c_action.read_input(mask, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

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

    c_action = c_action.Action_Center()
    c_action.read_input(command, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

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

    c_action = c_action.Action_RandomizeIons()
    c_action.read_input(command, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

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

    c_action = c_action.Action_ReplicateCell()
    c_action.read_input(command, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

    c_action.post_process()

    return mut_traj


def rotate_dihedral(traj=None, mask="", top=None):
    """rotate dihedral"""
    mut_traj = _assert_mutable(traj)

    c_action = c_action.Action_RotateDihedral()
    c_action.read_input(mask, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

    return mut_traj


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

    c_action = c_action.Action_SetDihedral()
    c_action.read_input(mask, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

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

    c_action = c_action.Action_SetVelocity()
    c_action.read_input(command, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

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

    c_action = c_action.Action_FixImagedBonds()
    c_action.read_input(mask, top=mut_traj.top)
    c_action.setup(mut_traj.top)

    for frame in mut_traj:
        c_action.compute(frame)

    return mut_traj


def atom_map(traj, ref, rmsfit=False):
    """compute atom mapping

    Parameters
    ----------
    traj : Trajectory-like
    ref : Trajectory-like
        reference trajectory
    rmsfit : bool, default False
        if True, perform RMS fitting

    Returns
    -------
    out : ndarray, mapping array

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> ref = pt.datafiles.load_tz2_ortho()
    >>> atom_indices = pt.atom_map(traj, ref)
    >>> assert len(atom_indices) == traj.top.n_atoms
    """
    if traj.n_atoms != ref.n_atoms:
        raise ValueError("number of atoms must be the same")

    # create mutable trajectory
    mut_traj = _assert_mutable(traj)
    mut_ref = _assert_mutable(ref) if hasattr(ref, 'xyz') else ref

    if rmsfit:
        fit_command = "rmsfit"
    else:
        fit_command = ""

    c_action = c_action.Action_AtomMap()
    c_action.read_input(fit_command, top=mut_traj.top)
    c_action.setup(mut_traj.top, crdinfo={'n_atoms': ref.top.n_atoms})

    # process reference
    for ref_frame in mut_ref:
        c_action.compute(ref_frame)

    # process target
    for frame in mut_traj:
        c_action.compute(frame)

    c_action.post_process()

    return c_action.atom_map


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
    c_action = c_action.Action_CheckChirality()
    c_action.read_input(mask, top=traj.top, dslist=c_dslist)
    c_action.setup(traj.top)

    for frame in traj:
        c_action.compute(frame)

    c_action.post_process()
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