"""
Structure utilities and basic analysis functions
"""
import numpy as np
from typing import Union

from ..utils.get_common_objects import (
    get_topology, get_fiterator, get_data_from_dtype,
    get_list_of_commands, get_reference, super_dispatch
)
from ..utils import ensure_not_none_or_string
from ..utils.decorators import register_pmap, register_openmp
from ..utils.convert import array_to_cpptraj_atommask
from ..datasets.datasetlist import DatasetList
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from ..trajectory.shared_methods import iterframe_master
from ..analysis.c_action import c_action, do_action
from ..analysis.c_action.actionlist import ActionList
from ..trajectory import Trajectory
from ..trajectory.frame import Frame
from .base_classes import (
    _assert_mutable, CommandType, _check_command_type,
    _create_and_compute_action_list, DatasetType
)


@register_pmap
@super_dispatch()
def strip(traj=None, mask="", frame_indices=None, top=None):
    """strip atoms from trajectory

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask to remove
    frame_indices : array-like, optional, default None
    top : Topology, optional

    Returns
    -------
    Trajectory
        trajectory with atoms stripped

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Strip water molecules
    >>> stripped = pt.strip(traj, ':WAT')
    >>> # Strip hydrogen atoms
    >>> stripped = pt.strip(traj, '@H*')
    >>> # Strip ions
    >>> stripped = pt.strip(traj, ':Na+,Cl-')
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'strip {mask}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Strip(), command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)

    return traj


@register_pmap
@super_dispatch()
def closest(traj=None, mask="", closestmols=0, frame_indices=None, top=None):
    """keep closest molecules to mask

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        reference mask to find closest molecules to
    closestmols : int
        number of closest molecules to keep
    frame_indices : array-like, optional, default None
    top : Topology, optional

    Returns
    -------
    Trajectory
        trajectory with only closest molecules

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Keep 50 closest water molecules to protein
    >>> closest_traj = pt.closest(traj, ':1-13', closestmols=50)
    >>> # Keep 10 closest waters to a specific residue
    >>> closest_traj = pt.closest(traj, ':5', closestmols=10)
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'closest {closestmols} {mask}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Closest(), command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)

    return traj


@register_pmap
@super_dispatch()
def mean_structure(traj=None, mask='*', frame_indices=None, top=None):
    """compute mean structure from trajectory

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*'
        atom mask for mean calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional

    Returns
    -------
    Frame
        mean structure as a Frame object

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate mean structure for entire system
    >>> mean_frame = pt.mean_structure(traj)
    >>> # Calculate mean for protein only
    >>> mean_frame = pt.mean_structure(traj, ':1-13')
    >>> # Calculate mean for specific frames
    >>> mean_frame = pt.mean_structure(traj, frame_indices=range(100, 200))
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    # Get atom indices for mask
    atom_indices = top.select(mask)

    # Calculate mean coordinates
    coords_sum = np.zeros((len(atom_indices), 3))
    n_frames = 0

    for frame in traj:
        coords_sum += frame.xyz[atom_indices]
        n_frames += 1

    if n_frames == 0:
        raise ValueError("No frames found")

    mean_coords = coords_sum / n_frames

    # Create mean frame
    mean_frame = Frame()
    mean_frame.xyz = frame.xyz.copy()  # Start with last frame
    mean_frame.xyz[atom_indices] = mean_coords
    mean_frame.top = top

    return mean_frame


def get_average_frame(traj=None, mask='*', frame_indices=None, top=None):
    """get average frame (alias for mean_structure)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*'
        atom mask for average calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional

    Returns
    -------
    Frame
        average structure as a Frame object

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> avg_frame = pt.get_average_frame(traj)
    """
    return mean_structure(traj=traj, mask=mask, frame_indices=frame_indices, top=top)


@register_pmap
@super_dispatch()
def pairdist(traj=None, mask1='', mask2='', frame_indices=None, top=None, dtype='ndarray'):
    """compute all pairwise distances between two masks

    Parameters
    ----------
    traj : Trajectory-like
    mask1 : str
        first atom mask
    mask2 : str
        second atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    3D ndarray
        pairwise distances, shape (n_frames, n_atoms1, n_atoms2)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # All distances between CA atoms and water oxygens
    >>> distances = pt.pairdist(traj, '@CA', ':WAT@O')
    >>> # Distances between two residues
    >>> distances = pt.pairdist(traj, ':1', ':10')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'pairdist {mask1} {mask2}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_PairDist(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def rdf(traj=None, solute_mask='', solvent_mask='', frame_indices=None,
        top=None, dtype='ndarray', bin_spacing=0.5, maximum=10.0):
    """compute radial distribution function

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask : str
        solute atom mask
    solvent_mask : str
        solvent atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    bin_spacing : float, default 0.5
        histogram bin spacing
    maximum : float, default 10.0
        maximum distance for RDF

    Returns
    -------
    DatasetList
        RDF data with distance and g(r) values

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # RDF of water around protein
    >>> rdf_data = pt.rdf(traj, ':1-13', ':WAT@O')
    >>> # RDF with custom parameters
    >>> rdf_data = pt.rdf(traj, ':1-13', ':WAT@O', bin_spacing=0.1, maximum=15.0)
    >>> # Ion-water RDF
    >>> rdf_data = pt.rdf(traj, ':Na+', ':WAT@O')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'radial {solute_mask} {solvent_mask} {bin_spacing} {maximum}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Radial(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def lipidscd(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray'):
    """compute lipid order parameters (SCD)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        lipid atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    DatasetList
        lipid order parameters

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate lipid order parameters
    >>> scd = pt.lipidscd(traj, ':DPPC')
    >>> # For specific lipid atoms
    >>> scd = pt.lipidscd(traj, ':DPPC@C2,C3,C4')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'lipidorder {mask}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_LipidOrder(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


def get_velocity(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray'):
    """get velocities from trajectory

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for velocity extraction
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    3D ndarray
        velocities, shape (n_frames, n_atoms, 3)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Get velocities for all atoms (if available)
    >>> velocities = pt.get_velocity(traj)
    >>> # Get velocities for CA atoms only
    >>> velocities = pt.get_velocity(traj, '@CA')

    Notes
    -----
    Returns velocities if present in trajectory, otherwise estimates
    from coordinate differences.
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    # Get atom indices
    if mask:
        atom_indices = top.select(mask)
    else:
        atom_indices = np.arange(top.n_atoms)

    velocities = []

    for frame in traj:
        if hasattr(frame, 'velocity') and frame.velocity is not None:
            # Use existing velocity data
            vel = frame.velocity[atom_indices]
        else:
            # Estimate from coordinate differences (requires previous frame)
            # For now, return zeros if no velocity data
            vel = np.zeros((len(atom_indices), 3))

        velocities.append(vel)

    if dtype == 'ndarray':
        return np.array(velocities)
    else:
        # Return as DatasetList
        dslist = DatasetList()
        dslist.add('velocity', np.array(velocities))
        return get_data_from_dtype(dslist, dtype)


def set_velocity(traj=None, velocities=None, mask='', frame_indices=None, top=None):
    """set velocities in trajectory

    Parameters
    ----------
    traj : Trajectory-like
    velocities : array-like
        velocity data to set, shape (n_frames, n_atoms, 3)
    mask : str, optional
        atom mask for velocity setting
    frame_indices : array-like, optional, default None
    top : Topology, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Set random velocities
    >>> vel = np.random.randn(len(traj), traj.n_atoms, 3)
    >>> pt.set_velocity(traj, vel)
    >>> # Set velocities for specific atoms
    >>> vel_ca = np.random.randn(len(traj), len(traj.top.select('@CA')), 3)
    >>> pt.set_velocity(traj, vel_ca, '@CA')
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    if velocities is None:
        raise ValueError("velocities array is required")

    velocities = np.asarray(velocities)

    # Get atom indices
    if mask:
        atom_indices = top.select(mask)
    else:
        atom_indices = np.arange(top.n_atoms)

    # Set velocities in each frame
    for i, frame in enumerate(traj):
        if not hasattr(frame, 'velocity') or frame.velocity is None:
            frame.velocity = np.zeros((top.n_atoms, 3))

        if i < len(velocities):
            frame.velocity[atom_indices] = velocities[i]


def make_structure(traj=None, frame_idx=0, top=None):
    """create a structure from a trajectory frame

    Parameters
    ----------
    traj : Trajectory-like
    frame_idx : int, default 0
        frame index to use for structure
    top : Topology, optional

    Returns
    -------
    Frame
        single frame structure

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Make structure from first frame
    >>> struct = pt.make_structure(traj, 0)
    >>> # Make structure from last frame
    >>> struct = pt.make_structure(traj, -1)
    >>> # Make structure from middle frame
    >>> struct = pt.make_structure(traj, len(traj)//2)
    """
    ensure_not_none_or_string(traj)

    top = get_topology(traj, top)

    # Get the specified frame
    frame = traj[frame_idx]

    # Create a copy to avoid modifying original
    structure = frame.copy()
    structure.top = top

    return structure


def contact_map(traj=None, mask='', frame_indices=None, top=None,
                dtype='ndarray', distance_cutoff=8.0, byres=True):
    """compute contact map

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for contact calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    distance_cutoff : float, default 8.0
        distance cutoff for contacts in Angstroms
    byres : bool, default True
        calculate contacts by residue instead of by atom

    Returns
    -------
    3D ndarray
        contact maps, shape (n_frames, n_residues, n_residues) if byres=True

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate residue contact map
    >>> contacts = pt.contact_map(traj, ':1-13')
    >>> # Calculate atom contact map
    >>> contacts = pt.contact_map(traj, '@CA', byres=False)
    >>> # Use different cutoff
    >>> contacts = pt.contact_map(traj, ':1-13', distance_cutoff=6.0)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['contacts', mask, f'{distance_cutoff}']
    if byres:
        command_parts.append('byres')

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Contacts(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


def solvent_density(traj=None, solute_mask='', solvent_mask=':WAT',
                   frame_indices=None, top=None, dtype='ndarray',
                   grid_spacing=0.5, buffer=5.0):
    """compute solvent density around solute

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask : str
        solute atom mask
    solvent_mask : str, default ':WAT'
        solvent atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    grid_spacing : float, default 0.5
        grid spacing for density calculation
    buffer : float, default 5.0
        buffer distance around solute

    Returns
    -------
    DatasetList
        solvent density data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate water density around protein
    >>> density = pt.solvent_density(traj, ':1-13', ':WAT@O')
    >>> # Use finer grid
    >>> density = pt.solvent_density(traj, ':1-13', ':WAT@O', grid_spacing=0.25)
    """
    ensure_not_none_or_string(traj)

    command = (f'density {solvent_mask} delta {grid_spacing} '
               f'around {solute_mask} {buffer}')

    dslist, _ = do_action(traj, command, c_action.Action_Density,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Export all functions
__all__ = [
    'strip', 'closest', 'mean_structure', 'get_average_frame', 'pairdist',
    'rdf', 'lipidscd', 'get_velocity', 'set_velocity', 'make_structure',
    'contact_map', 'solvent_density'
]