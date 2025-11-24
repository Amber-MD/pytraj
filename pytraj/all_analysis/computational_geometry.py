"""
Computational geometry and algorithm tools
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
from .base_classes import (
    CommandType, _check_command_type, _create_and_compute_action_list,
    DatasetType
)


@register_pmap
@super_dispatch()
def search_neighbors(traj=None, mask='', cutoff=5.0, frame_indices=None,
                    top=None, dtype='dataset'):
    """search for neighboring atoms within cutoff distance

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for neighbor search
    cutoff : float, default 5.0
        distance cutoff for neighbors
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'

    Returns
    -------
    DatasetList
        neighbor lists and distances

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Find neighbors within 5.0 Å of CA atoms
    >>> neighbors = pt.search_neighbors(traj, '@CA', cutoff=5.0)
    >>> # Find neighbors within 3.0 Å of a specific residue
    >>> neighbors = pt.search_neighbors(traj, ':5', cutoff=3.0)
    >>> # Find water neighbors around protein
    >>> neighbors = pt.search_neighbors(traj, ':1-13', cutoff=4.0)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'neighbors {mask} {cutoff}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Neighbors(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def hausdorff(traj=None, ref=None, mask='', frame_indices=None, top=None,
              dtype='ndarray'):
    """calculate Hausdorff distance between structures

    Parameters
    ----------
    traj : Trajectory-like
    ref : Frame or int
        reference structure
    mask : str, optional
        atom mask for Hausdorff calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    1D ndarray
        Hausdorff distances for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Hausdorff distance to first frame
    >>> hausdorff_dist = pt.hausdorff(traj, ref=0, mask='@CA')
    >>> # Hausdorff distance for all atoms
    >>> hausdorff_dist = pt.hausdorff(traj, ref=traj[0])
    >>> # Hausdorff distance for specific residues
    >>> hausdorff_dist = pt.hausdorff(traj, ref=0, mask=':1-10')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    # Get reference structure
    if ref is None:
        ref_frame = traj[0]
    elif isinstance(ref, int):
        ref_frame = traj[ref]
    else:
        ref_frame = ref

    command_parts = ['hausdorff', 'reference']
    if mask:
        command_parts.extend([mask, 'refmask', mask])

    command = ' '.join(command_parts)

    # Setup reference
    dslist = CpptrajDatasetList()
    dslist.add(DatasetType.REFERENCE, name='hausdorff_ref')
    dslist[0].top = ref_frame.top or top
    dslist[0].add_frame(ref_frame)

    action_list = ActionList()
    action_list.add(c_action.Action_Hausdorff(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def atom_map(traj=None, ref=None, mask='', frame_indices=None, top=None,
             dtype='dataset', algorithm='hungarian'):
    """perform atom mapping between structures

    Parameters
    ----------
    traj : Trajectory-like
    ref : Frame or int
        reference structure for mapping
    mask : str, optional
        atom mask for mapping
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    algorithm : str, default 'hungarian'
        mapping algorithm ('hungarian', 'greedy')

    Returns
    -------
    DatasetList
        atom mapping results

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Map atoms to first frame
    >>> mapping = pt.atom_map(traj, ref=0, mask='@CA')
    >>> # Use greedy algorithm
    >>> mapping = pt.atom_map(traj, ref=0, mask='@CA', algorithm='greedy')
    >>> # Map specific residues
    >>> mapping = pt.atom_map(traj, ref=0, mask=':1-5@CA,CB,CG')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    # Get reference structure
    if ref is None:
        ref_frame = traj[0]
    elif isinstance(ref, int):
        ref_frame = traj[ref]
    else:
        ref_frame = ref

    command_parts = ['atommap', 'reference', algorithm]
    if mask:
        command_parts.extend([mask, 'refmask', mask])

    command = ' '.join(command_parts)

    # Setup reference
    dslist = CpptrajDatasetList()
    dslist.add(DatasetType.REFERENCE, name='atommap_ref')
    dslist[0].top = ref_frame.top or top
    dslist[0].add_frame(ref_frame)

    action_list = ActionList()
    action_list.add(c_action.Action_AtomMap(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def lowestcurve(traj=None, mask='', frame_indices=None, top=None,
                dtype='ndarray', npoints=100):
    """calculate lowest curve through 3D points

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for curve calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    npoints : int, default 100
        number of points for curve discretization

    Returns
    -------
    DatasetList
        lowest curve data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate lowest curve through CA atoms
    >>> curve = pt.lowestcurve(traj, '@CA')
    >>> # Use more points for smoother curve
    >>> curve = pt.lowestcurve(traj, '@CA', npoints=200)
    >>> # Calculate for backbone atoms
    >>> curve = pt.lowestcurve(traj, '@N,CA,C')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'lowestcurve {mask} npoints {npoints}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_LowestCurve(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def crank(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray'):
    """perform crankshaft move analysis

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for crankshaft analysis
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    DatasetList
        crankshaft analysis results

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Crankshaft analysis on backbone
    >>> crank_data = pt.crank(traj, '@N,CA,C')
    >>> # Analysis on specific region
    >>> crank_data = pt.crank(traj, ':5-10@N,CA,C')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'crank {mask}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Crank(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


def convex_hull(traj=None, mask='', frame_indices=None, top=None, dtype='dataset'):
    """calculate convex hull of atomic coordinates

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for convex hull calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'

    Returns
    -------
    DatasetList
        convex hull vertices and properties

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Convex hull of CA atoms
    >>> hull = pt.convex_hull(traj, '@CA')
    >>> # Convex hull of entire protein
    >>> hull = pt.convex_hull(traj, ':1-13')
    """
    ensure_not_none_or_string(traj)

    command = f'hull {mask}'

    dslist, _ = do_action(traj, command, c_action.Action_Hull,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def closest_distance_between_sets(traj=None, mask1='', mask2='',
                                 frame_indices=None, top=None, dtype='ndarray'):
    """find closest distance between two atom sets

    Parameters
    ----------
    traj : Trajectory-like
    mask1 : str
        first atom set mask
    mask2 : str
        second atom set mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    1D ndarray
        closest distances for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Closest distance between two residues
    >>> min_dist = pt.closest_distance_between_sets(traj, ':1', ':10')
    >>> # Closest distance between protein and ligand
    >>> min_dist = pt.closest_distance_between_sets(traj, ':1-13', ':LIG')
    """
    ensure_not_none_or_string(traj)

    command = f'mindist {mask1} {mask2}'

    dslist, _ = do_action(traj, command, c_action.Action_MinDist,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def spatial_clustering(traj=None, mask='', frame_indices=None, top=None,
                      dtype='dataset', algorithm='dbscan', eps=3.0, min_samples=3):
    """perform spatial clustering of atomic coordinates

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for clustering
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    algorithm : str, default 'dbscan'
        clustering algorithm
    eps : float, default 3.0
        maximum distance for clustering
    min_samples : int, default 3
        minimum samples per cluster

    Returns
    -------
    DatasetList
        clustering results

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Cluster water molecules
    >>> clusters = pt.spatial_clustering(traj, ':WAT@O')
    >>> # Use different parameters
    >>> clusters = pt.spatial_clustering(traj, ':WAT@O', eps=4.0, min_samples=5)
    """
    ensure_not_none_or_string(traj)

    command = f'cluster {mask} {algorithm} eps {eps} minsamples {min_samples}'

    dslist, _ = do_action(traj, command, c_action.Action_Cluster,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def distance_geometry_constraints(traj=None, mask='', frame_indices=None,
                                 top=None, dtype='dataset', lower_bound=2.0,
                                 upper_bound=8.0):
    """generate distance geometry constraints

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for constraint generation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    lower_bound : float, default 2.0
        lower distance bound
    upper_bound : float, default 8.0
        upper distance bound

    Returns
    -------
    DatasetList
        distance constraints

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Generate constraints for CA atoms
    >>> constraints = pt.distance_geometry_constraints(traj, '@CA')
    >>> # Custom bounds
    >>> constraints = pt.distance_geometry_constraints(traj, '@CA',
    ...                                               lower_bound=1.5, upper_bound=10.0)
    """
    ensure_not_none_or_string(traj)

    command = f'distgeom {mask} lower {lower_bound} upper {upper_bound}'

    dslist, _ = do_action(traj, command, c_action.Action_DistGeom,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Export all functions
__all__ = [
    'search_neighbors', 'hausdorff', 'atom_map', 'lowestcurve', 'crank',
    'convex_hull', 'closest_distance_between_sets', 'spatial_clustering',
    'distance_geometry_constraints'
]