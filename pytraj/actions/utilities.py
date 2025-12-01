"""
Utility and miscellaneous functions  
"""
from .base import *

__all__ = [
    'get_velocity', 'mean_structure', 'get_average_frame', 'multidihedral', 
    'search_neighbors', 'native_contacts', 'grid', 'transform', 'lowestcurve',
    'superpose', 'rotdif', 'lipidscd', 'xtalsymm', 'analyze_modes', 'ti',
    'hausdorff', 'permute_dihedrals'
]


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
        c_action = c_action.Action_Dihedral()
        c_action.read_input(commands[0], top=traj.top, dslist=c_dslist)
        c_action.setup(traj.top)

        for frame in traj:
            c_action.compute(frame)

        c_action.post_process()
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
    c_action = c_action.Action_SearchNeighbors()
    c_action.read_input(command, top=traj.top, dslist=c_dslist)
    c_action.setup(traj.top)

    for frame in traj:
        c_action.compute(frame)

    c_action.post_process()
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

    c_action = c_action.Action_NativeContacts()
    c_action.read_input(command + " ref ref", top=traj.top, dslist=c_dslist)
    c_action.setup(traj.top)

    for frame in traj:
        c_action.compute(frame)

    c_action.post_process()
    
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
    c_action = c_action.Action_Grid()
    c_action.read_input(command, top=traj.top, dslist=c_dslist)
    c_action.setup(traj.top)

    for frame in traj:
        c_action.compute(frame)

    c_action.post_process()
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
        c_action = c_action.Action_Transform()
        c_action.read_input(command, top=mut_traj.top)
        c_action.setup(mut_traj.top)

        for frame in mut_traj:
            c_action.compute(frame)
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
    c_analysis.Analysis_RotDif(command + " matrices", dslist=c_dslist)
    
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
    c_action = c_action.Action_LipidOrder()
    c_action.read_input(command, top=traj.top, dslist=c_dslist)
    c_action.setup(traj.top)

    for frame in traj:
        c_action.compute(frame)

    c_action.post_process()
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

    c_action = c_action.Action_XtalSymm()
    c_action.read_input(command, top=traj.top, dslist=c_dslist)
    c_action.setup(traj.top)

    for frame in traj:
        c_action.compute(frame)

    c_action.post_process()
    
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

    c_action = c_action.Action_PermuteDihedrals()
    c_action.read_input(command, top=mut_traj.top, dslist=c_dslist)
    c_action.setup(mut_traj.top)

    # permute
    c_action.compute(None)
    c_action.post_process()

    c_dslist._pop(0)
    c_dslist._pop(0)
    
    return mut_traj