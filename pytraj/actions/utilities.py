"""
Utility and miscellaneous functions
"""
from .base import *
from .base import _assert_mutable, add_reference_dataset
from ..trajectory.trajectory_iterator import TrajectoryIterator

__all__ = [
    'get_velocity', 'mean_structure', 'get_average_frame', 'multidihedral',
    'search_neighbors', 'native_contacts', 'grid', 'transform', 'lowestcurve',
    'superpose', 'rotdif', 'lipidscd', 'xtalsymm', 'analyze_modes', 'ti',
    'hausdorff', 'permute_dihedrals', 'in_voxel', 'count_in_voxel',
    'pair_distance', 'closest_atom', 'crank', 'check_structure', 'set_velocity'
]

# _assert_mutable is imported from base.py


# voxel center and xyz are tuples
def in_voxel(voxel_cntr, xyz, delta):
    return (xyz[0] >= voxel_cntr[0] - delta and xyz[0] <= voxel_cntr[0] + delta
           ) and (xyz[1] >= voxel_cntr[1] - delta and xyz[1] <= voxel_cntr[1] +
                  delta) and (xyz[2] >= voxel_cntr[2] - delta and
                              xyz[2] <= voxel_cntr[2] + delta)


@super_dispatch()
def count_in_voxel(traj=None, mask="", voxel_cntr=(0, 0, 0), voxel_size=5, frame_indices=None, top=None):
    """For a voxel with center xyz and size voxel_size, find atoms that match a given mask
    that are contained in that voxel over the course of a trajectory.

    This analysis command is meant to be used with GIST analysis to estimate water residence time
    for a given region. When running GIST on a trajectory, we get a list of voxels and their
    associated rotational/translational entropy. To analyze how long solvent molecules are
    residing in various regions of the macromolecular surface, we can compute residence times
    using the survival time correlation function. This can be done by plotting the number of
    solvent molecules that reside in the voxel at frame 1 and then plotting the decay in the
    number of those original solvent molecules as they diffuse out of the region and are
    replaced. This can provide insight into the behavior of solvent atoms closely interfacing
    with the protein or nucleic acid of study. The voxel_cntr is the center of the voxel with
    format [x,y,z] where x,y,z are in real space unit. Similarly, voxel_size should be in the
    same real space unit as the trajectory coordinates. The starting frame should be specified
    on the trajectory (often frame 1) and then analysis should calculate the residence time over
    remaining frames. This will result in the final output being 1 less value than the number
    of frames in the trajectory.

    Parameters
    ----------
    traj: Trajectory-like
    mask: string containing atoms to be searched for in the voxel. Default empty string is all atoms.
    voxel_cntr: list, default (0,0,0)
        center of voxel
    voxel_size: int, default 5
        size of voxel measured from center to edge.
    frame_indices : array-like, optional
        frame indices
    top : Topology, optional

    Returns
    -------
    1D ndarray, one dimensional numpy array containing the counts for each frame.
    """
    lives_in_voxel = []
    topology = get_topology(traj, top)
    population = topology.atom_indices(mask)
    delta = voxel_size / 2

    for frame in traj:
        frame_voxAtoms = []
        for atm in population:
            coord = frame.atom(atm)
            if (in_voxel(voxel_cntr, coord, delta)):
                frame_voxAtoms.append(atm)
        lives_in_voxel.append(frame_voxAtoms)

    return lives_in_voxel


def pair_distance(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)


@register_pmap
def closest_atom(top=None, frame=None, point=(0, 0, 0), mask=""):
    """for a given xyz coordinate in a frame, find the closest atom

    Parameters
    ----------
    top: Topology object
    frame: Frame object
    point: Tuple containing 3 elements, for x y and z coordinate of point.
    Default point is the origin
    mask: string containing atoms to find, in atom mask syntax. Default is empty string
    which contains all the atoms

    Returns
    -------
    Index of atom closest to point in xyz coordinate space.

    Notes
    -----
    Topology needs to contain atoms that match the atom mask passed in, and frame needs to have
    xyz coordinates for all atoms. Point should be a tuple of length 3 with format (x, y, z)

    Examples
    --------
    >>> import pytraj as pt
    >>> # find closest atom to origin in a given trajectory
    >>> traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
    >>> frame = traj[0]
    >>> pt.closest_atom(traj.top, traj[0], (0, 0, 0))
    205
    """

    if (len(top.atom_indices(mask)) == 0):
        raise ValueError(
            "Please pass in a topology file with atoms that match the mask in it"
        )

    closest_dist = None
    closest_idx = None
    atoms = top.atom_indices(mask)
    for atm in atoms:
        coord = frame.atom(atm)
        if ((closest_dist is None) or
            (pair_distance(coord, point) < closest_dist)):
            closest_dist = pair_distance(coord, point)
            closest_idx = atm

    return closest_idx


@register_pmap
def mean_structure(traj,
                   mask='',
                   frame_indices=None,
                   dtype='frame',
                   autoimage=False,
                   rmsfit=None,
                   top=None):
    '''get mean structure for a given mask and given frame_indices

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask : None or str, default None (all atoms)
    frame_indices : iterable that produces integer, default None, optional
        frame indices
    dtype: str, {'frame', 'traj'}, default 'frame'
        return type, either Frame (does not have Topology information) or 'traj'
    autoimage : bool, default False
        if True, performa autoimage
    rmsfit : object, {Frame, int, tuple, None}, default None
        if rmsfit is not None, perform rms fit to reference.

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # get average structure from whole traj, all atoms
    >>> frame = pt.mean_structure(traj)

    >>> # get average structure from several frames, '@CA' atoms"
    >>> frame = pt.mean_structure(traj, '@CA', frame_indices=range(2, 8, 2))

    >>> # get average structure but do autoimage and rmsfit to 1st Frame
    >>> frame = pt.mean_structure(traj(autoimage=True, rmsfit=0))

    >>> # get average structure but do autoimage and rmsfit to 1st Frame.
    >>> frame = pt.mean_structure(traj(autoimage=True, rmsfit=0, frame_indices=[0, 5, 6]))

    Notes
    -----
    if autoimage=True and having rmsfit, perform autoimage first and do rmsfit
    '''
    # note: we not yet use @super_dispatch due to extra 'rmsfit'
    # TODO: do it.
    topology = get_topology(traj, top)
    try:
        frame_iterator = traj.iterframe(autoimage=autoimage,
                                        rmsfit=rmsfit,
                                        frame_indices=frame_indices)
    except AttributeError:
        frame_iterator = get_fiterator(traj, frame_indices)

    action_datasets = CpptrajDatasetList()
    if not isinstance(mask, str):
        mask = array_to_cpptraj_atommask(mask)

    # add "crdset s1" to trick cpptraj dump coords to DatSetList
    command = mask + " crdset s1"

    action_average = c_action.Action_Average()
    action_average(command, frame_iterator, topology, dslist=action_datasets)

    # need to call this method so cpptraj will write
    action_average.post_process()

    frame = action_datasets[0].get_frame()
    if dtype.lower() == 'frame':
        return frame
    elif dtype.lower() in ['traj', 'trajectory']:
        new_topology = topology if mask == '' else topology[mask]
        return Trajectory(xyz=frame.xyz.reshape(1, frame.n_atoms, 3).copy(),
                          top=new_topology)
    else:
        raise ValueError('dtype must be frame or trajectory')


# create alias
get_average_frame = mean_structure


def get_velocity(traj, mask=None, top=None, frame_indices=None):
    '''get velocity for specify frames with given mask

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask : str, default None (use all atoms), optional
        atom mask
    top : Topology, optional
    frame_indices : iterable that produces integer, default None, optional
        if not None, only get velocity for given frame indices

    Returns
    -------
    out : 3D numpy array, shape (n_frames, n_atoms, 3)

    Examples
    --------
    >>> vels = pt.get_velocity(traj, frame_indices=[0, 3]) # doctest: +SKIP

    Notes
    -----
    Since pytraj has limited support for force and velocity info, if user wants to
    load both from disk, should iterate the TrajectoryIterator (got by pytraj.iterload method)

    .. code-block:: python

        import pytraj as pt
        forces = []
        velocities = []

        traj = pt.iterload(filename, topology_filename)

        for frame in traj:
            forces.append(frame.force)
            velocities.append(frame.velocity)

        # Note: pytraj can efficiently load arbitary frame indices to memory
        for frame in pt.iterframe(traj, frame_indices=[0, 8, 8, 100, 1000]): pass
    '''
    if mask is None:
        atm_indices = None
    else:
        if not isinstance(mask, str):
            # array-like
            atm_indices = mask
        else:
            atm_indices = traj.top.select(mask)

    fi = traj.iterframe(frame_indices=frame_indices)
    n_atoms = traj.n_atoms if mask is None else len(atm_indices)
    n_frames = fi.n_frames

    data = np.empty((n_frames, n_atoms, 3), dtype='f8')
    for idx, frame in enumerate(fi):
        if not frame.has_velocity():
            raise ValueError('frame does not have velocity')
        data[idx] = frame.velocity if mask is None else frame.velocity[
            atm_indices]
    return data


@super_dispatch()
def multidihedral(traj=None,
                  dihedral_types=None,
                  resrange=None,
                  define_new_type=None,
                  range360=False,
                  dtype='dataset',
                  top=None,
                  frame_indices=None,
                  mass=False):
    """perform dihedral search

    Parameters
    ----------
    traj : Trajectory-like object
    dihedral_types : dihedral type, default None
        if None, calculate all supported dihedrals
    resrange : str | array-like
        residue range for searching. If `resrange` is string, use index starting with 1
        (cpptraj convertion)
        if `resrange` is array-like, use index starting from 0 (python convention)
    define_new_type : str
        define new type for searching
    range360 : bool, default False
        if True: use 0-360
    top : Topology | str, optional
        only need to have 'top' if can not find it in `traj`


    Returns
    -------
    pytraj.DatasetList (use `values` attribute to get raw `numpy` array)

    Notes
    -----
        Dataset lables show residue number in 1-based index

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.multidihedral(traj)
    >>> data = pt.multidihedral(traj, 'phi psi')
    >>> data = pt.multidihedral(traj, resrange=range(8))
    >>> data = pt.multidihedral(traj, range360=True)
    >>> data = pt.multidihedral(traj, resrange='1,3,5')
    >>> data = pt.multidihedral(traj, dihedral_types='phi psi')
    >>> data = pt.multidihedral(traj, dihedral_types='phi psi', resrange='3-7')
    >>> data = pt.multidihedral(traj, dihedral_types='phi psi', resrange=[3, 4, 8])
    """
    # Process resrange
    resrange_str = None
    if resrange:
        if isinstance(resrange, str):
            resrange_str = str(resrange)
        else:
            from pytraj.utils import convert as cv
            resrange_str = str(cv.array_to_cpptraj_range(resrange))

    # Process define_new_type
    define_new_type_str = None
    if define_new_type:
        define_new_type_str = f"dihtype {str(define_new_type)}"

    # Build command using CommandBuilder
    command = (CommandBuilder().add(
        str(dihedral_types), condition=dihedral_types is not None).add(
            "resrange", resrange_str, condition=resrange_str
            is not None).add(define_new_type_str,
                             condition=define_new_type_str
                             is not None).add("range360",
                                              condition=range360).add("mass",
                                              condition=mass).build())

    action_datasets, _ = do_action(traj, command, c_action.Action_MultiDihedral)
    return get_data_from_dtype(action_datasets, dtype=dtype)


@super_dispatch()
def search_neighbors(traj=None,
                     mask='',
                     distance=3.0,
                     frame_indices=None,
                     dtype='dataset',
                     top=None):
    """find neighbors within a cutoff

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    distance : float, default 3.0
        cutoff distance
    frame_indices : array-like, optional
    dtype : str, default 'dataset'
        return data type
    top : Topology, optional

    Returns
    -------
    out : DatasetList
        neighbors information
    """
    action_datasets = DatasetList()
    topology = get_topology(traj, top)

    # Handle frame_indices by filtering frames
    if frame_indices is not None:
        selected_frames = []
        for idx, frame in enumerate(iterframe_master(traj)):
            if idx in frame_indices:
                selected_frames.append((idx, frame))
        frame_list = selected_frames
    else:
        frame_list = [
            (idx, frame) for idx, frame in enumerate(iterframe_master(traj))
        ]

    for frame_idx, frame in frame_list:
        topology.set_reference(frame)
        selected_indices = topology.select(mask)
        action_datasets.append({str(frame_idx): np.asarray(selected_indices)})

    return get_data_from_dtype(action_datasets, dtype)


@super_dispatch(refindex=3)
def native_contacts(traj=None,
                    mask="",
                    mask2="",
                    ref=0,
                    dtype='dataset',
                    distance=7.0,
                    mindist=None,
                    maxdist=None,
                    image=True,
                    include_solvent=False,
                    byres=False,
                    series=False,
                    first=False,
                    frame_indices=None,
                    options='',
                    top=None):
    """Compute native contacts.

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default ""
        First atom selection
    mask2 : str, default ""
        Second atom selection
    ref : {Frame, int}, default 0
        Reference frame or index
    dtype : str, default 'dataset'
        Return data type
    distance : float, default 7.0
        Distance cutoff for contacts (Angstroms)
    mindist : float, optional
        Minimum distance cutoff
    maxdist : float, optional
        Maximum distance cutoff (overrides distance if set)
    image : bool, default True
        Apply periodic boundary conditions
    include_solvent : bool, default False
        Include solvent atoms in contact calculation
    byres : bool, default False
        Calculate contacts per residue pair
    series : bool, default False
        Calculate contact time series for each pair
    first : bool, default False
        Use first occurrence for contact identification
    frame_indices : array-like, optional
        Frame indices to analyze
    options : str
        Extra cpptraj command(s).
    top : Topology, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # use 1st frame as reference, don't need specify ref Frame
    >>> data = pt.native_contacts(traj)

    >>> # explicitly specify reference, specify distance cutoff
    >>> ref = traj[3]
    >>> data = pt.native_contacts(traj, ref=ref, distance=8.0)

    >>> # use integer array for mask
    >>> data = pt.native_contacts(traj, mask=range(100), mask2=[200, 201], ref=ref, distance=8.0)

    >>> # Enhanced parameters: per-residue contacts with distance range
    >>> data = pt.native_contacts(traj, ref=ref, byres=True, mindist=3.0, maxdist=10.0)

    >>> # Time series analysis with first contact identification
    >>> data = pt.native_contacts(traj, ref=ref, series=True, first=True)
    """
    ref = get_reference(traj, ref)
    native_contacts_action = c_action.Action_NativeContacts()
    action_datasets = CpptrajDatasetList()

    if not isinstance(mask2, str):
        # [1, 3, 5] to "@1,3,5
        mask2 = array_to_cpptraj_atommask(mask2)
    mask_str = ' '.join((mask, mask2))

    # Build command with enhanced parameters
    cmd_builder = (CommandBuilder().add("ref myframe").add(mask_str))

    # Distance parameters - maxdist overrides distance if set
    if maxdist is not None:
        cmd_builder.add("maxdist", str(maxdist))
    else:
        cmd_builder.add("distance", str(distance))

    if mindist is not None:
        cmd_builder.add("mindist", str(mindist))

    # Boolean options
    cmd_builder.add("noimage", condition=not image)
    cmd_builder.add("includesolvent", condition=include_solvent)
    cmd_builder.add("byresidue", condition=byres)
    cmd_builder.add("series", condition=series)
    cmd_builder.add("first", condition=first)
    cmd_builder.add(options, condition=bool(options))

    command = cmd_builder.build()

    add_reference_dataset(action_datasets, 'myframe', ref, top)
    native_contacts_action(command, traj, top=top, dslist=action_datasets)
    action_datasets.remove_at(0)

    return get_data_from_dtype(action_datasets, dtype=dtype)


@super_dispatch()
def grid(traj=None, command="", top=None, dtype='dataset', frame_indices=None):
    """perform grid analysis

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        cpptraj grid command
    top : Topology, optional
    dtype : str, default 'dataset'
        return data type
    frame_indices : array-like, optional
        frame indices

    Returns
    -------
    out : DatasetList
    """
    action_datasets, _ = do_action(traj, command, c_action.Action_Grid, top=top)
    return get_data_from_dtype(action_datasets, dtype=dtype)


def transform(traj, by, top=None, frame_indices=None):
    '''transform pytraj.Trajectory by a series of cpptraj's commands

    Parameters
    ----------
    traj : Mutable Trajectory
    by : list of cpptraj commands
    top : Topology, optional
    frame_indices : {None, array-like}, default None
        if not None, perform tranformation for specific frames.

    Returns
    -------
    transformed Trajectory. Trajectory's coordinates will be inplace-updated

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # perform 'autoimage', then 'rms', then 'center'
    >>> traj = pt.transform(traj[:], by=['autoimage', 'rms', 'center :1-5'])
    '''
    return traj.transform(by, frame_indices=frame_indices)


def lowestcurve(data, points=10, step=0.2):
    '''compute lowest curve for data

    Parameters
    ----------
    data : 2D array-like
    points : number of lowest points in each bin, default 10
    step : step size, default 0.2

    Returns
    -------
    2d array
    '''
    command = f'mydata points {points} step {step}'

    data = np.asarray(data).T

    runner = AnalysisRunner(c_analysis.Analysis_LowestCurve)
    runner.add_dataset(DatasetType.XYMESH, 'mydata', data)

    runner.run_analysis(command)

    return np.array(
        [runner.datasets[-1]._xcrd(),
         np.array(runner.datasets[-1].values)])


def superpose(traj, *args, **kwd):
    """alias of rmsd.superpose method"""
    from .rmsd import superpose as _superpose
    return _superpose(traj, *args, **kwd)


def rotdif(matrices, command):
    """perform rotational diffusion analysis

    Parameters
    ----------
    matrices : array-like
        rotation matrices
    command : str
        cpptraj rotdif command

    Returns
    -------
    out : DatasetList
    """
    matrices = np.asarray(matrices, dtype='f8')

    c_dslist = CpptrajDatasetList()

    # add matrix data
    matrix_dataset = c_dslist.add('matrix3x3', 'matrices')
    for mat in matrices:
        matrix_dataset._append_from_array(mat.flatten())

    # run analysis
    c_analysis.Analysis_Rotdif(command + " matrices", dslist=c_dslist)

    return c_dslist


@super_dispatch()
def lipidscd(traj, mask='', options='', dtype='dict', top=None, frame_indices=None):
    """compute lipid order parameters

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    options : str, optional
        extra options
    dtype : str, default 'dict'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional
        frame indices

    Returns
    -------
    out : dict or DatasetList
    """
    command = mask + " " + options
    action_datasets, _ = do_action(traj, command, c_action.Action_LipidOrder, top=top)
    return get_data_from_dtype(action_datasets, dtype=dtype)


@super_dispatch()
def xtalsymm(traj, mask='', options='', ref=None, frame_indices=None, **kwargs):
    """compute crystal symmetry analysis

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    options : str, optional
        extra options
    ref : Frame, optional
        reference frame
    frame_indices : array-like, optional
        frame indices
    **kwargs : additional options

    Returns
    -------
    out : DatasetList
    """
    command = mask + " " + options
    for key, val in kwargs.items():
        command += f" {key} {val}"

    c_dslist = CpptrajDatasetList()

    if ref is not None:
        ref_frame = get_reference(traj, ref)
        add_reference_dataset(c_dslist, 'ref', ref_frame, ref_frame.top or traj.top)
        command += " ref ref"

    do_action(traj, command, c_action.Action_XtalSymm, dslist=c_dslist)

    if ref is not None:
        c_dslist.remove_at(0)  # remove reference

    return c_dslist


def analyze_modes(mode_type,
                  eigenvectors,
                  eigenvalues,
                  scalar_type='mwcovar',
                  options='',
                  dtype='dict'):
    runner = AnalysisRunner(c_analysis.Analysis_Modes)
    my_modes = 'my_modes'
    runner.add_dataset(DatasetType.MODES, my_modes, None)

    modes = runner.datasets[-1]
    modes.scalar_type = scalar_type
    modes._allocate_avgcoords(eigenvectors.shape[1])
    modes._set_modes(False, eigenvectors.shape[0], eigenvectors.shape[1],
                     eigenvalues, eigenvectors.flatten())

    command = ' '.join((mode_type, 'name {}'.format(my_modes), options))
    runner.run_analysis(command)

    runner.datasets.remove_at(0)
    return get_data_from_dtype(runner.datasets, dtype=dtype)


def ti(fn, options=''):
    """compute thermodynamic integration

    Parameters
    ----------
    fn : str
        filename
    options : str, optional
        extra options

    Returns
    -------
    out : dict
        TI results
    """
    command = f"readdata {fn} {options}"

    c_dslist = CpptrajDatasetList()
    c_analysis.Analysis_TI(command, dslist=c_dslist)

    return get_data_from_dtype(c_dslist, dtype='dict')


def hausdorff(matrix, options='', dtype='ndarray'):
    """
    Parameters
    ----------
    matrix : 2D array
    option : str
        cpptraj options

    Returns
    -------
    out : 1D numpy array

    Notes
    -----
        - cpptraj help: pytraj.info('hausdorff')
    """
    runner = AnalysisRunner(c_analysis.Analysis_Hausdorff)
    runner.add_dataset(DatasetType.MATRIX_DBL, "my_matrix", matrix)
    runner.run_analysis(f"my_matrix {options}")
    runner.datasets.remove_at(0)  # Remove input matrix
    return get_data_from_dtype(runner.datasets, dtype)


def permute_dihedrals(traj, filename, options='', top=None):
    """
    Parameters
    ----------
    traj : Trajectory like
    filename : str
        Output filename for resulted trajectory
    options: str
        cpptraj's option. Do not specify `outtraj` here since
        it's specified in `filename`.

    This function returns None.
    """
    from ..core.c_core import CpptrajState, Command
    from .base import DatasetType

    state = CpptrajState()

    top_data = state.data.add(DatasetType.TOPOLOGY, name='my_top')
    top_data.data = traj.top

    ref_data = state.data.add(DatasetType.COORDS, name='my_coords')
    ref_data.top = traj.top
    for frame in traj:
        ref_data.add_frame(frame)

    command = f'permutedihedrals crdset my_coords {options} outtraj {filename}'

    with Command() as executor:
        executor.dispatch(state, command)

    state.data.remove_at(0)
    state.data.remove_at(0)


def check_structure(traj,
                    mask='',
                    options='',
                    frame_indices=None,
                    top=None,
                    dtype='ndarray'):
    """check if the structure is ok or not

    Parameters
    ----------
    traj : Trajectory-like
    mask: str, default all atoms
    options : str, default ''
        extra cpptraj options
    dtype : str, default 'ndarray'

    Returns
    -------
    out : Tuple[np.ndarray, str]
        number of problems for each frame and detail

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_rna()
    >>> failures = pt.check_structure(traj[:1])
    """
    command = ' '.join((mask, options))
    action_datasets, c_stdout = do_action(traj, command,
                                          c_action.Action_CheckStructure)
    return get_data_from_dtype(action_datasets, dtype=dtype), c_stdout


def crank(data0, data1, mode='distance', dtype='ndarray'):
    """Crank-shaft analysis

    Parameters
    ----------
    data0 : array-like
    data1 : array-like
    mode : str, {'distance', 'angle'}
    dtype : str

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> distances = pt.distance(traj, [':3 :7', ':8 :12'])
    >>> out = pt.crank(distances[0], distances[1])

    Notes
    -----
    Same as `crank` in cpptraj
    """
    runner = AnalysisRunner(c_analysis.Analysis_CrankShaft)
    runner.add_dataset(DatasetType.DOUBLE, "d0", data0)
    runner.add_dataset(DatasetType.DOUBLE, "d1", data1)

    command = ' '.join((mode, 'd0', 'd1'))
    with capture_stdout() as (out, err):
        runner.run_analysis(command)
    return out.read()


def set_velocity(traj, temperature=298, ig=10, options='', top=None):
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
