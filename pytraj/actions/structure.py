"""
Structure analysis functions: RMSD, surfaces, volumes, pucker analysis
"""
from .base import *

__all__ = [
    'radgyr', 'radgyr_tensor', 'surf', 'molsurf', 'volume', 'pucker', 'multipucker',
    'watershell', 'rmsf', 'bfactors', '_calc_vector_center', 'center_of_mass', 'center_of_geometry'
]


def _calc_vector_center(traj=None,
                        mask="",
                        top=None,
                        mass=False,
                        dtype='ndarray',
                        frame_indices=None):

    action_datasets = CpptrajDatasetList()
    action_datasets.set_own_memory(False)  # need this to avoid segmentation fault
    execute_vector_action = c_action.Action_Vector()
    command = "center " + mask

    if mass:
        command += " mass"

    # Ensure top is a valid Topology object
    if top is None:
        top = traj.top

    execute_vector_action(command, traj, top=top, dslist=action_datasets)
    return get_data_from_dtype(action_datasets, dtype=dtype)


def center_of_mass(traj=None,
                   mask='',
                   top=None,
                   dtype='ndarray',
                   frame_indices=None):
    '''compute center of mass

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # compute center of mass residue 3 for first 2 frames.
    array([[-0.661702  ,  6.69124347,  3.35159413],
           [ 0.5620708 ,  7.82263042, -0.72707798]])
    '''
    # note: do not use super_dispatch for this method since
    # we already use for _calc_vector_center
    return _calc_vector_center(
        traj=traj,
        mask=mask,
        top=top,
        mass=True,
        dtype=dtype,
        frame_indices=frame_indices)


@register_pmap
@super_dispatch()
def center_of_geometry(traj=None,
                       mask="",
                       top=None,
                       dtype='ndarray',
                       frame_indices=None):

    atom_mask_obj = top(mask)
    action_datasets = CpptrajDatasetList()
    action_datasets.add(DatasetType.VECTOR)

    for frame in iterframe_master(traj):
        action_datasets[0].append(frame.center_of_geometry(atom_mask_obj))
    return get_data_from_dtype(action_datasets, dtype=dtype)


@super_dispatch()
def radgyr(traj=None,
           mask='',
           dtype='ndarray',
           top=None,
           frame_indices=None):
    """compute radgyr

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape (n_frames,)
    """
    action = c_action.Action_RadGyr()
    c_dslist = CpptrajDatasetList()
    action.read_input(mask, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def radgyr_tensor(traj=None, mask='', dtype='ndarray', top=None, frame_indices=None):
    """compute radius of gyration tensor

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray
    """
    action = c_action.Action_RadGyr()
    c_dslist = CpptrajDatasetList()
    tensor_command = mask + " tensor"
    action.read_input(tensor_command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def surf(traj=None, mask="", dtype='ndarray', frame_indices=None, top=None):
    """compute solvent accessible surface area

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape (n_frames,)
    """
    action = c_action.Action_Surf()
    c_dslist = CpptrajDatasetList()
    action.read_input(mask, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def molsurf(traj=None,
            mask='',
            probe=1.4,
            radii_type='',
            dtype='ndarray',
            frame_indices=None,
            top=None):
    """compute molecular surface area

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    probe : float, default 1.4
        probe radius
    radii_type : str, optional
        VDW radii to use for calculation
    dtype : str, default 'ndarray'
        return data type
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    out : ndarray, shape (n_frames,)
    """
    command = f"{mask} probe {probe}"
    if radii_type:
        command += f" {radii_type}"

    action = c_action.Action_Molsurf()
    c_dslist = CpptrajDatasetList()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def volume(traj=None, mask="", top=None, dtype='ndarray', frame_indices=None):
    """compute volume

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
    top : Topology, optional
    dtype : str, default 'ndarray'
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape (n_frames,)
    """
    action = c_action.Action_Volume()
    c_dslist = CpptrajDatasetList()
    action.read_input(mask, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def watershell(traj=None,
               mask='',
               solvent_mask=':WAT@O',
               lower=3.4,
               upper=5.0,
               dtype='ndarray',
               frame_indices=None,
               top=None):
    """compute waterhell

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for finding nearest atoms
    solvent_mask : str, default ':WAT@O'
        solvent mask
    lower : float, default 3.4
        lower bound distance
    upper : float, default 5.0
        upper bound distance
    dtype : str, default 'ndarray'
        return data type
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    out : ndarray, shape (n_frames,)
    """
    command = f"{solvent_mask} around {mask} lower {lower} upper {upper}"
    action = c_action.Action_Watershell()
    c_dslist = CpptrajDatasetList()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def rmsf(traj=None,
         mask="",
         top=None,
         dtype='ndarray',
         frame_indices=None,
         options=''):
    '''compute atomicfluct (RMSF)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str or 1D-array
        atom mask. If not given, use all atoms
    options : str, additional cpptraj options ('byres', 'bymask', 'byatom', 'calcadp')

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.rmsf(traj, '@CA') # or pt.atomicfluct
    >>> data[:3]
    array([[  5.        ,   0.61822273],
           [ 16.        ,   0.5627449 ],
           [ 40.        ,   0.53717119]])
    '''
    command = ' '.join((mask, options))
    action_datasets, _ = do_action(traj, command, c_action.Action_AtomicFluct)
    return get_data_from_dtype(action_datasets, dtype=dtype)


def bfactors(traj=None,
             mask="",
             byres=True,
             top=None,
             dtype='ndarray',
             frame_indices=None):
    """calculate pseudo bfactor

    Notes
    -----
    This is **NOT** getting bfactor from xray, but computing bfactor from simulation.

    Parameters
    ----------
    traj: Trajectory-like
    mask: str, mask

    Returns
    -------
    if dtype is 'ndarray' (default), return a numpy array
    with shape=(n_atoms/n_residues, 2) ([atom_or_residue_idx, value])

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2')
    >>> traj = pt.load(fn, tn, mask='!:WAT')
    >>> traj = pt.superpose(traj)
    >>> bfactor = pt.bfactors(traj, byres=True)
    """
    byres_text = "byres" if byres else ""

    # need to convert to string mask
    # do not use super_dispatch again
    if not isinstance(mask, str):
        mask = array_to_cpptraj_atommask(mask)
    command_ = " ".join((mask, byres_text, "bfactor"))
    return rmsf(
        traj=traj,
        mask=command_,
        top=top,
        dtype=dtype,
        frame_indices=frame_indices)


def pucker(traj=None,
           mask=None,
           method="altona",
           dtype='ndarray',
           top=None,
           frame_indices=None):
    """compute pucker

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    method : str, default 'altona'
        pucker calculation method. Either 'altona' or 'cremer'
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape (n_frames,)
    """
    command = mask or ""
    if method == 'cremer':
        command += " cremer"
    elif method == 'altona':
        command += " altona"
    else:
        raise ValueError("method must be either 'altona' or 'cremer'")

    action = c_action.Action_Pucker()
    c_dslist = CpptrajDatasetList()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def multipucker(traj=None, resrange=None, method="altona", range360=False, amplitude=False, ampout=None, theta=False, thetaout=None, offset=None, out=None, puckertype=None, extra_options="", dtype='dict', top=None, frame_indices=None):
    """Perform multi-pucker analysis.

    Parameters
    ----------
    traj : Trajectory-like
    resrange : str or array-like, optional
        Residue range for the analysis. If not provided, all residues are used.
    method : str, default "altona"
        Pucker calculation method. Options: "altona" or "cremer".
    range360 : bool, default False
        Use 0-360 degree range if True, otherwise use -180 to 180.
    amplitude : bool, default False
        Include amplitude in the output if True.
    ampout : str, optional
        Output filename for amplitude values.
    theta : bool, default False
        Include theta values in the output if True.
    thetaout : str, optional
        Output filename for theta values.
    offset : float, optional
        Offset to add to pucker values (in degrees).
    out : str, optional
        Output filename for the results.
    puckertype : str, optional
        Specific pucker type to calculate (e.g., "furanoid:C2:C3:C4:C5:O2", "pyranoid:C1:C2:C3:C4:C5:O5", "pyranose")
    extra_options : str, optional
        Additional cpptraj options for future extensibility.
    dtype : str, default 'dict'
        Output data type.
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    dict or DatasetList depending on dtype
    """
    if resrange:
        if isinstance(resrange, str):
            resrange_str = resrange
        else:
            from ..utils.convert import array_to_cpptraj_range
            resrange_str = array_to_cpptraj_range(resrange)
    else:
        resrange_str = None

    command = (CommandBuilder()
               .add("resrange", resrange_str, condition=resrange_str is not None)
               .add("puckertype", puckertype, condition=puckertype is not None)
               .add(method)
               .add("range360", condition=range360)
               .add("amplitude", condition=amplitude)
               .add("ampout", ampout, condition=ampout is not None)
               .add("theta", condition=theta)
               .add("thetaout", thetaout, condition=thetaout is not None)
               .add("offset", str(offset), condition=offset is not None)
               .add("out", out, condition=out is not None)
               .add(extra_options, condition=bool(extra_options))
               .build())

    action_datasets, _ = do_action(traj, command, c_action.Action_MultiPucker)
    return get_data_from_dtype(action_datasets, dtype=dtype)