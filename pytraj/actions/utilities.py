"""
Utility and miscellaneous functions
"""
from .base import *
from .base import _assert_mutable
from ..trajectory.trajectory_iterator import TrajectoryIterator


__all__ = [
    'get_velocity', 'mean_structure', 'get_average_frame', 'multidihedral',
    'search_neighbors', 'native_contacts', 'grid', 'transform', 'lowestcurve',
    'superpose', 'rotdif', 'lipidscd', 'xtalsymm', 'analyze_modes', 'ti',
    'hausdorff', 'permute_dihedrals', 'in_voxel',
    'count_in_voxel', 'pair_distance', 'closest_atom', 'crank', 'check_structure',
    'set_velocity'
]

# _assert_mutable is imported from base.py


# voxel center and xyz are tuples
def in_voxel(voxel_cntr, xyz, delta):
    return (xyz[0] >= voxel_cntr[0] - delta and xyz[0] <= voxel_cntr[0] +
            delta) and (xyz[1] >= voxel_cntr[1] - delta
                        and xyz[1] <= voxel_cntr[1] + delta) and (
                            xyz[2] >= voxel_cntr[2] - delta
                            and xyz[2] <= voxel_cntr[2] + delta)


@register_pmap
def count_in_voxel(traj=None, mask="", voxel_cntr=(0, 0, 0), voxel_size=5):
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

    Returns
    -------
    1D ndarray, one dimensional numpy array containing the counts for each frame.
    """
    lives_in_voxel = []
    population = traj.top.atom_indices(mask)
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
    >>> traj = pt.load("test.nc", "test.parm7")
    >>> atom_index = pt.closest_atom(traj.top, traj[0], (1.0, 2.0, 0.5))
    """
    indices = top.atom_indices(mask)
    min_distance = float('inf')
    closest_index = -1

    for i in indices:
        x = frame.xyz[i][0]
        y = frame.xyz[i][1]
        z = frame.xyz[i][2]

        distance = pair_distance(point, (x, y, z))
        if distance < min_distance:
            min_distance = distance
            closest_index = i

    return closest_index


def mean_structure(traj,
                   mask='*',
                   frame_indices=None,
                   dtype='trajectory',
                   top=None):
    """compute mean structure

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*'
        atom mask
    frame_indices : array-like, optional
        frame indices
    dtype : str, default 'trajectory'
        return data type
    top : Topology, optional

    Returns
    -------
    out : Trajectory (with 1 frame) or Frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> meanframe = pt.mean_structure(traj)
    >>> meanframe.n_frames
    1
    """
    if frame_indices is not None:
        traj = traj[frame_indices]

    if mask == '*':
        xyz = np.mean(traj.xyz, axis=0)
    else:
        atom_indices = traj.top.select(mask)
        xyz = np.mean(traj.xyz[:, atom_indices], axis=0)

    if dtype == 'trajectory':
        # create new trajectory with mean coordinates
        mean_traj = traj[:1].copy()
        if mask == '*':
            mean_traj.xyz[0] = xyz
        else:
            mean_traj.xyz[0, atom_indices] = xyz
        return mean_traj
    elif dtype == 'frame':
        mean_frame = traj[0].copy()
        if mask == '*':
            mean_frame.xyz = xyz
        else:
            mean_frame.xyz[atom_indices] = xyz
        return mean_frame
    else:
        return xyz


# create alias
get_average_frame = mean_structure


def get_velocity(traj, mask=None, frame_indices=None):
    """get velocity

    Parameters
    ----------
    traj : Trajectory-like
    mask : str or array-like, optional
        if given, use this mask to select atoms
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape (n_frames, n_atoms, 3)

    Notes
    -----
    This method will return 3D array with shape=(n_frames, n_atoms, 3).
    If there is no velocity information, return None
    """
    if hasattr(traj, 'xyz_v'):
        if mask is None:
            return traj.xyz_v
        else:
            if isinstance(mask, str):
                atom_indices = traj.top.select(mask)
            else:
                atom_indices = mask
            return traj.xyz_v[:, atom_indices]
    else:
        return None


@super_dispatch()
def multidihedral(traj=None,
                  mask='',
                  resrange=None,
                  dihedral_types=None,
                  dtype='dataset',
                  top=None,
                  frame_indices=None):
    """compute multiple dihedrals

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    resrange : array-like, optional
        residue range
    dihedral_types : array-like, optional
        dihedral types to calculate
    dtype : str, default 'dataset'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : DatasetList or ndarray
    """
    if dihedral_types is None:
        dihedral_types = ['phi', 'psi']
    if resrange is None:
        resrange = range(1, traj.top.n_residues + 1)

    commands = []
    for resid in resrange:
        for dih_type in dihedral_types:
            commands.append(f"dihedral :{resid}@{dih_type}")

    if len(commands) == 1:
        # single dihedral
        c_dslist = CpptrajDatasetList()
        action = c_action.Action_Dihedral()
        action.read_input(commands[0], top=traj.top, dslist=c_dslist)
        action.setup(traj.top)

        for frame in traj:
            action.compute(frame)

        action.post_process()
        return get_data_from_dtype(c_dslist, dtype=dtype)
    else:
        # multiple dihedrals
        from .geometry import _create_and_compute_action_list
        return _create_and_compute_action_list(commands, traj,
                                             c_action.Action_Dihedral, dtype=dtype)


@super_dispatch()
def search_neighbors(traj=None,
                     mask='',
                     distance=3.0,
                     frame_indices=None,
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
    top : Topology, optional

    Returns
    -------
    out : dict
        neighbors information
    """
    command = f"{mask} {distance}"

    c_dslist = CpptrajDatasetList()
    action = c_action.Action_SearchNeighbors()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype='dict')


@super_dispatch(refindex=3)
def native_contacts(traj=None,
                    mask='',
                    mask2='',
                    ref=None,
                    distance=7.0,
                    frame_indices=None,
                    top=None):
    """compute native contacts

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        first mask
    mask2 : str, optional
        second mask
    ref : int or Frame, default None
        reference for native contacts, default first frame
    distance : float, default 7.0
        cutoff distance
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    out : ndarray, shape (n_frames,)
        fraction of native contacts
    """
    if ref is None:
        ref = traj[0]

    command = f"{mask} {mask2} distance {distance}"

    # setup reference
    c_dslist = CpptrajDatasetList()

    # add reference frame
    ref_frame = get_reference(traj, ref)
    ref_dataset = c_dslist.add('reference', 'ref')
    ref_dataset.top = ref_frame.top or traj.top
    ref_dataset.add_frame(ref_frame)

    action = c_action.Action_NativeContacts()
    action.read_input(command + " ref ref", top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()

    # remove reference dataset
    c_dslist._pop(0)
    return get_data_from_dtype(c_dslist, dtype='ndarray')


@super_dispatch()
def grid(traj=None, command="", top=None, dtype='dataset'):
    """perform grid analysis

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        cpptraj grid command
    top : Topology, optional
    dtype : str, default 'dataset'
        return data type

    Returns
    -------
    out : DatasetList
    """
    c_dslist = CpptrajDatasetList()
    action = c_action.Action_Grid()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


def transform(traj, by, frame_indices=None):
    """transform trajectory coordinates

    Parameters
    ----------
    traj : Trajectory-like
    by : str or array-like
        transformation matrix or cpptraj transform command
    frame_indices : array-like, optional

    Returns
    -------
    traj : transformed trajectory
    """
    mut_traj = _assert_mutable(traj)

    if isinstance(by, str):
        command = by
        action = c_action.Action_Transform()
        action.read_input(command, top=mut_traj.top)
        action.setup(mut_traj.top)

        for frame in mut_traj:
            action.compute(frame)
    else:
        # assume matrix transformation
        by = np.asarray(by, dtype='f8')
        if by.shape != (4, 4):
            raise ValueError("transformation matrix must be 4x4")

        for frame in mut_traj:
            # apply transformation matrix to coordinates
            xyz_homo = np.column_stack([frame.xyz, np.ones(frame.n_atoms)])
            frame.xyz = xyz_homo.dot(by.T)[:, :3]

    return mut_traj


def lowestcurve(data, points=10, step=0.2):
    """find lowest curve

    Parameters
    ----------
    data : array-like
        input data
    points : int, default 10
        number of points
    step : float, default 0.2
        step size

    Returns
    -------
    out : ndarray
        lowest curve data
    """
    data = np.asarray(data, dtype='f8')

    c_dslist = CpptrajDatasetList()
    dataset = c_dslist.add('double', 'data')
    dataset.data = data

    command = f"lowestcurve data points {points} step {step}"
    c_analysis.Analysis_LowestCurve(command, dslist=c_dslist)

    return c_dslist[-1].values


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
def lipidscd(traj, mask='', options='', dtype='dict'):
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

    Returns
    -------
    out : dict or DatasetList
    """
    command = mask + " " + options

    c_dslist = CpptrajDatasetList()
    action = c_action.Action_LipidOrder()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def xtalsymm(traj, mask='', options='', ref=None, **kwargs):
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
        ref_dataset = c_dslist.add('reference', 'ref')
        ref_dataset.top = ref_frame.top or traj.top
        ref_dataset.add_frame(ref_frame)
        command += " ref ref"

    action = c_action.Action_XtalSymm()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()

    if ref is not None:
        c_dslist._pop(0)  # remove reference

    return c_dslist


def analyze_modes(mode_type,
                  eigenvalues,
                  eigenvectors,
                  scalar_type='mwcovar',
                  options='',
                  mask='*'):
    """analyze normal modes

    Parameters
    ----------
    mode_type : str
        type of analysis (fluct, displ, corr)
    eigenvalues : array-like
        eigenvalues
    eigenvectors : array-like
        eigenvectors
    scalar_type : str, default 'mwcovar'
        scalar type
    options : str, optional
        extra options
    mask : str, default '*'
        atom mask

    Returns
    -------
    out : DatasetList
    """
    eigenvalues = np.asarray(eigenvalues, dtype='f8')
    eigenvectors = np.asarray(eigenvectors, dtype='f8')

    c_dslist = CpptrajDatasetList()

    # add eigenvalues
    eval_dataset = c_dslist.add('modes', 'eigenvalues')
    eval_dataset.data = eigenvalues

    # add eigenvectors
    evec_dataset = c_dslist.add('modes', 'eigenvectors')
    evec_dataset.data = eigenvectors

    command = f"{mode_type} eigenvalues eigenvectors {scalar_type} {mask} {options}"

    # run analysis
    c_analysis.Analysis_Modes(command, dslist=c_dslist)

    return c_dslist


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
    """compute Hausdorff distance

    Parameters
    ----------
    matrix : array-like
        distance matrix
    options : str, optional
        extra options
    dtype : str, default 'ndarray'
        return data type

    Returns
    -------
    out : ndarray or DatasetList
    """
    matrix = np.asarray(matrix, dtype='f8')

    c_dslist = CpptrajDatasetList()

    # add matrix
    matrix_dataset = c_dslist.add('matrix_dbl', 'matrix')
    matrix_dataset.data = matrix

    command = f"hausdorff matrix {options}"

    # run analysis
    c_analysis.Analysis_Hausdorff(command, dslist=c_dslist)

    return get_data_from_dtype(c_dslist, dtype=dtype)


def permute_dihedrals(traj, filename, options=''):
    """permute dihedral angles

    Parameters
    ----------
    traj : Trajectory-like
    filename : str
        output filename
    options : str, optional
        extra options

    Returns
    -------
    traj : Trajectory with permuted dihedrals
    """
    mut_traj = _assert_mutable(traj)

    command = f"crdout {filename} {options}"

    # create data to hold trajectory
    c_dslist = CpptrajDatasetList()
    coords_data = c_dslist.add('coords', name='_DEFAULTCRD_')
    coords_data.top = mut_traj.top
    for frame in mut_traj:
        coords_data.append(frame)

    action = c_action.Action_PermuteDihedrals()
    action.read_input(command, top=mut_traj.top, dslist=c_dslist)
    action.setup(mut_traj.top)

    # permute
    action.compute(None)
    action.post_process()

    c_dslist._pop(0)
    c_dslist._pop(0)

    return mut_traj


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


@super_dispatch()
def set_velocity(traj, temperature=298, ig=10, options=''):
    """set velocity for atoms

    Parameters
    ----------
    traj: Trajectory object
    temperature: float, default 298 K
        temperature for Maxwell-Boltzmann distribution
    ig: int, default 10
        random seed
    options: str, optional
        additional options

    Returns
    -------
    None (modify `traj` inplace)
    """
    if hasattr(traj, 'top'):
        command = f"setvelocity temp {temperature} ig {ig} {options}"
        do_action(traj, command, c_action.Action_SetVelocity)
