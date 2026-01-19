"""
Structure analysis functions: RMSD, surfaces, volumes, pucker analysis
"""
from .base import *

__all__ = [
    'radgyr', 'radgyr_tensor', 'surf', 'molsurf', 'volume', 'pucker',
    'multipucker', 'watershell', 'rmsf', 'bfactors', '_calc_vector_center',
    'center_of_mass', 'center_of_geometry'
]


def _calc_vector_center(traj=None,
                        mask="",
                        top=None,
                        mass=False,
                        dtype='ndarray',
                        frame_indices=None):

    action_datasets = CpptrajDatasetList()
    action_datasets.set_own_memory(
        False)  # need this to avoid segmentation fault
    execute_vector_action = c_action.Action_Vector()
    command = "center " + mask

    if mass:
        command += " mass"

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
    return _calc_vector_center(traj=traj,
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


@register_pmap
@super_dispatch()
def radgyr(traj=None,
           mask="",
           top=None,
           nomax=True,
           mass=True,
           tensor=False,
           frame_indices=None,
           dtype='ndarray'):
    '''compute radius of gyration

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default ""
        atom selection
    top : Topology, optional
    nomax : bool, default True
        do not calculate tensor eigenvalues maximum
    mass : bool, default True
        use mass-weighted calculation
    tensor : bool, default False
        calculate radius of gyration tensor
    frame_indices : array-like, optional
    dtype : str, default 'ndarray'
        return data type

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.radgyr(traj, '@CA')
    >>> data = pt.radgyr(traj, '!:WAT', nomax=False)
    >>> data = pt.radgyr(traj, '@CA', frame_indices=[2, 4, 6])
    >>> # Mass-weighted with tensor calculation
    >>> data = pt.radgyr(traj, '@CA', mass=True, tensor=True)
    '''
    command_parts = [mask]
    if nomax:
        command_parts.append('nomax')
    if not mass:
        command_parts.append('geom')  # Use geometric center instead of mass-weighted
    if tensor:
        command_parts.append('tensor')
    command = " ".join(command_parts)
    action_datasets, _ = do_action(traj, command, c_action.Action_Radgyr)
    return get_data_from_dtype(action_datasets, dtype)


@super_dispatch()
def radgyr_tensor(traj=None,
                  mask="",
                  top=None,
                  frame_indices=None,
                  dtype='ndarray'):
    '''compute radius of gyration with tensor

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default ""
        atom selection
    top : Topology, optional
    frame_indices : array-like, optional
    dtype : str, default 'ndarray'
        return data type

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.radgyr_tensor(traj, '@CA',)
    >>> data = pt.radgyr_tensor(traj, '@CA', frame_indices=[2, 4, 6])

    Returns
    -------
    out : Dict[str, np.ndarray]
    '''
    nomax_ = 'nomax'
    command = " ".join((mask, nomax_, "tensor"))
    action_datasets, _ = do_action(traj, command, c_action.Action_Radgyr)
    k0, v0 = action_datasets[0].key, action_datasets[0].values.copy(
    )  # use copy to avoid early memory free
    k1, v1 = action_datasets[1].key, action_datasets[1].possible_data6
    if dtype == 'dict':
        return {k0: v0, k1: v1}
    elif dtype == 'ndarray':
        # Return structured ndarray or tuple for backward compatibility
        # Note: returning tuple for backward compatibility, could be changed in future version
        return v0, v1
    else:
        return get_data_from_dtype(action_datasets, dtype=dtype)


@super_dispatch()
def surf(traj=None, mask="", dtype='ndarray', frame_indices=None, top=None,
         offset=1.4, nbrcut=2.5, solutemask=None):
    """calc surf (LCPO method) - compute solvent accessible surface area

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    dtype : str, default 'ndarray'
        return data type
    frame_indices : array-like, optional
    top : Topology, optional
    offset : float, default 1.4
        van der Waals offset in Angstroms
    nbrcut : float, default 2.5
        Cutoff for determining neighbors in Angstroms
    solutemask : str, optional
        Mask to define solute atoms

    Returns
    -------
    out : ndarray, shape (n_frames,)
        surface area values for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.surf(traj, '@CA')
    """
    if traj is None:
        raise ValueError('trajectory is required')

    # Build command with parameters
    command = mask
    if offset != 1.4:
        command += f" offset {offset}"
    if nbrcut != 2.5:
        command += f" nbrcut {nbrcut}"
    if solutemask is not None:
        command += f" solutemask {solutemask}"

    action_datasets, _ = do_action(traj, command, c_action.Action_Surf)
    return get_data_from_dtype(action_datasets, dtype=dtype)


@super_dispatch()
def molsurf(traj=None,
            mask="",
            probe=1.4,
            offset=0.0,
            dtype='ndarray',
            frame_indices=None,
            top=None):
    '''calc molsurf

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> pt.molsurf(traj, '@CA')
    array([ 458.51409637,  459.64784573,  456.54690793,  467.72939574,
            462.45908781,  458.70327554,  454.40514806,  455.15015576,
            468.70566447,  456.0058624 ])
    >>> pt.molsurf(traj, '!:WAT')
    array([ 1079.1395679 ,  1090.79759341,  1069.65127413,  1096.0810919 ,
            1091.65862234,  1091.68906298,  1085.53105392,  1069.22510187,
            1079.70803583,  1075.8151414 ])
    '''
    probe_value = 'probe ' + str(probe)
    offset_value = 'offset ' + str(offset) if offset != 0. else ''
    command = ' '.join((mask, probe_value, offset_value))
    action_datasets, _ = do_action(traj, command, c_action.Action_Molsurf)
    return get_data_from_dtype(action_datasets, dtype)


@super_dispatch()
def volume(traj=None, mask="", top=None, dtype='ndarray', frame_indices=None):
    """compute volume

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
    dtype : str, default 'ndarray'
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    out : ndarray, shape (n_frames,)
        volume values for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> vol = pt.volume(traj, '@CA')
    """
    if traj is None:
        raise ValueError('trajectory is required')

    action_datasets, _ = do_action(traj, mask, c_action.Action_Volume)
    return get_data_from_dtype(action_datasets, dtype=dtype)


@register_pmap
@register_openmp
@super_dispatch()
def watershell(traj=None,
               solute_mask='',
               solvent_mask=':WAT',
               lower=3.4,
               upper=5.0,
               image=True,
               dtype='dataset',
               frame_indices=None,
               top=None):
    """(adapted from cpptraj doc): Calculate numbers of waters in 1st and 2nd solvation shells
    (defined by <lower cut> (default 3.4 Ang.) and <upper cut> (default 5.0 Ang.)

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask: solute mask
    solvent_mask: solvent mask
    lower : double, default 3.4
        lower cut distance
    upper : double, default 5.0
        upper cut distance
    image : bool, defaul True
        do autoimage if True
    dtype : return type, defaul 'dataset'
    top : Topology, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.watershell(traj, solute_mask='!:WAT')
    >>> data = pt.watershell(traj, solute_mask='!:WAT', lower=5.0, upper=10.)
    """
    if solute_mask in [None, '']:
        raise ValueError('must provide solute mask')

    command = (CommandBuilder().add(solute_mask).add("lower", str(lower)).add(
        "upper", str(upper)).add("noimage",
                                 condition=not image).add(solvent_mask,
                                                          condition=solvent_mask
                                                          is not None).build())

    action_datasets, _ = do_action(traj, command, c_action.Action_Watershell)
    return get_data_from_dtype(action_datasets, dtype=dtype)


@super_dispatch()
def rmsf(traj=None,
         mask="",
         top=None,
         dtype='ndarray',
         frame_indices=None,
         byres=False,
         byatom=True,
         bymask=False,
         calcadp=False,
         options=''):
    '''compute atomicfluct (RMSF)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str or 1D-array
        atom mask. If not given, use all atoms
    byres : bool, default False
        Calculate fluctuations per residue
    byatom : bool, default True
        Calculate fluctuations per atom (default behavior)
    bymask : bool, default False
        Calculate single value for entire mask
    calcadp : bool, default False
        Calculate atomic displacement parameters
    options : str, additional cpptraj options

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.rmsf(traj, '@CA') # or pt.atomicfluct
    >>> data[:3]
    array([[  5.        ,   0.61822273],
           [ 16.        ,   0.5627449 ],
           [ 40.        ,   0.53717119]])

    >>> # Enhanced parameters example
    >>> # Calculate RMSF per residue
    >>> data_byres = pt.rmsf(traj, ':1-10', byres=True)
    >>> # Calculate single RMSF value for entire mask
    >>> data_bymask = pt.rmsf(traj, '@CA', bymask=True)
    '''
    # Build command with granularity options
    command_parts = [mask]

    # Granularity options - mutually exclusive
    if byres:
        command_parts.append('byres')
    elif bymask:
        command_parts.append('bymask')
    # byatom is default, no need to add flag

    # Additional options
    if calcadp:
        command_parts.append('calcadp')
    if options:
        command_parts.append(options)

    command = ' '.join(command_parts)
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
    return rmsf(traj=traj,
                mask=command_,
                top=top,
                dtype=dtype,
                frame_indices=frame_indices)


def pucker(traj=None,
           pucker_mask=("C1'", "C2'", "C3'", "C4'", "O4'"),
           resrange=None,
           top=None,
           dtype='dataset',
           range360=False,
           method='altona',
           use_com=True,
           amplitude=False,
           offset=None):
    """compute pucker

    Parameters
    ----------
    traj : Trajectory-like
    pucker_mask : str
    resrange : None or array of int
    top : Topology, optional
    dtype : str, return type
    range360: bool, use 360 or 180 scale
    method : {'altona', 'cremer'}, default 'altona'
    use_com : bool
    amplitude : bool, default False
    offset : None or float

    Returns
    -------
    Dataset
    """
    top_ = get_topology(traj, top)
    if resrange is None:
        resrange = range(top_.n_residues)

    c_dslist = CpptrajDatasetList()

    for res in resrange:
        atom_mask = " ".join(
            (":" + str(res + 1) + '@' + x for x in pucker_mask))
        name = "pucker_res" + str(res + 1)

        command = (CommandBuilder().add(name).add(atom_mask).add(
            "range360", condition=range360).add(method).add(
                "geom", condition=not use_com).add(
                    "amplitude", condition=amplitude).add("offset",
                                                          str(offset),
                                                          condition=offset
                                                          is not None).build())

        act = c_action.Action_Pucker()
        act(command, traj, top=top_, dslist=c_dslist)

    return get_data_from_dtype(c_dslist, dtype)


@super_dispatch()
def multipucker(traj=None,
                resrange=None,
                method="altona",
                range360=False,
                amplitude=False,
                ampout=None,
                theta=False,
                thetaout=None,
                offset=None,
                out=None,
                puckertype=None,
                extra_options="",
                dtype='dict',
                top=None,
                frame_indices=None):
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

    command = (CommandBuilder().add(
        "resrange", resrange_str, condition=resrange_str is not None).add(
            "puckertype", puckertype, condition=puckertype
            is not None).add(method).add("range360", condition=range360).add(
                "amplitude", condition=amplitude).add(
                    "ampout", ampout, condition=ampout
                    is not None).add("theta", condition=theta).add(
                        "thetaout", thetaout, condition=thetaout
                        is not None).add(
                            "offset", str(offset), condition=offset
                            is not None).add(
                                "out", out, condition=out is not None).add(
                                    extra_options,
                                    condition=bool(extra_options)).build())

    action_datasets, _ = do_action(traj, command, c_action.Action_MultiPucker)
    return get_data_from_dtype(action_datasets, dtype=dtype)
