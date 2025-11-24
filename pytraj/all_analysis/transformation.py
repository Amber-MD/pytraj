"""
Coordinate transformation functions
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
from .base_classes import (
    _assert_mutable, CommandType, _check_command_type,
    _create_and_compute_action_list, DatasetType
)


@register_pmap
@super_dispatch()
def translate(traj=None, command="", top=None, frame_indices=None):
    """translate coordinate

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        translation command (e.g., 'x 5.0' or 'center :1')
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Translate by 5.0 in x direction
    >>> pt.translate(traj, 'x 5.0')
    >>> # Translate to center a residue at origin
    >>> pt.translate(traj, 'center :1')
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    full_command = f'translate {command}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Translate(), full_command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def rotate(traj=None, command="", frame_indices=None, top=None):
    """rotate coordinates

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        rotation command (e.g., 'x 90' for 90 degrees around x-axis)
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Rotate 90 degrees around x-axis
    >>> pt.rotate(traj, 'x 90')
    >>> # Rotate around arbitrary axis
    >>> pt.rotate(traj, 'x 45 y 30 z 60')
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    full_command = f'rotate {command}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Rotate(), full_command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def center(traj=None, mask="", center='box', mass=False,
           top=None, frame_indices=None):
    """Center coordinates in mask to specified point

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask to center
    center : str, default 'box'
        center type ('box', 'origin', or coordinate values)
    mass : bool, default False
        if True, use mass-weighted centering
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Center protein in box
    >>> pt.center(traj, ':1-13', center='box')
    >>> # Center at origin
    >>> pt.center(traj, '@CA', center='origin')
    >>> # Mass-weighted centering
    >>> pt.center(traj, ':1-13', mass=True)
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['center', mask]
    if mass:
        command_parts.append('mass')
    command_parts.append(center)

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Center(), command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def scale(traj=None, command="", frame_indices=None, top=None):
    """scale coordinates

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        scaling command (e.g., 'x 1.5' or 'xyz 2.0')
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Scale x coordinates by 1.5
    >>> pt.scale(traj, 'x 1.5')
    >>> # Scale all coordinates by 2.0
    >>> pt.scale(traj, 'xyz 2.0')
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    full_command = f'scale {command}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Scale(), full_command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def align(traj=None, mask='', ref=None, ref_mask='', mass=False,
          frame_indices=None, top=None):
    """align trajectory to reference

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for alignment
    ref : Frame or int, default None
        reference frame
    ref_mask : str, optional
        reference mask (if different from mask)
    mass : bool, default False
        if True, use mass-weighted alignment
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Align backbone to first frame
    >>> pt.align(traj, '@CA,C,N')
    >>> # Align to specific reference
    >>> pt.align(traj, '@CA', ref=traj[0])
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)
    ref_mask = ref_mask or mask

    # Get reference
    if ref is None:
        ref_frame = traj[0]
    elif hasattr(ref, 'xyz'):
        ref_frame = ref
    else:
        ref_frame = traj[ref]

    # Build command
    command_parts = ['rms', 'reference', mask]
    if ref_mask != mask:
        command_parts.extend(['refmask', ref_mask])
    if mass:
        command_parts.append('mass')

    command = ' '.join(command_parts)

    # Setup datasets and reference
    dslist = CpptrajDatasetList()
    dslist.add(DatasetType.REFERENCE, name='myref')
    dslist[0].top = ref_frame.top or top
    dslist[0].add_frame(ref_frame)

    action_list = ActionList()
    action_list.add(c_action.Action_Rmsd(), command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def superpose(traj=None, mask='', ref=None, ref_mask='', mass=False,
              frame_indices=None, top=None):
    """superpose trajectory to reference (alias for align)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for superposition
    ref : Frame or int, default None
        reference frame
    ref_mask : str, optional
        reference mask (if different from mask)
    mass : bool, default False
        if True, use mass-weighted superposition
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> pt.superpose(traj, '@CA')
    """
    return align(traj=traj, mask=mask, ref=ref, ref_mask=ref_mask,
                 mass=mass, frame_indices=frame_indices, top=top)


@register_pmap
@super_dispatch()
def autoimage(traj=None, command="", frame_indices=None, top=None):
    """automatically image coordinates to primary unit cell

    Parameters
    ----------
    traj : Trajectory-like
    command : str, optional
        additional autoimage options
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> pt.autoimage(traj)
    >>> # With options
    >>> pt.autoimage(traj, 'anchor :1')
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    full_command = f'autoimage {command}'.strip()

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_AutoImage(), full_command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def image(traj=None, command="", frame_indices=None, top=None):
    """image coordinates

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        image command options
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Image with center
    >>> pt.image(traj, 'center :1-13')
    >>> # Image specific residues
    >>> pt.image(traj, 'byres :WAT')
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    full_command = f'image {command}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Image(), full_command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def transform(traj=None, command="", frame_indices=None, top=None):
    """apply coordinate transformation

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        transformation matrix or command
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Apply transformation matrix
    >>> pt.transform(traj, 'matrix rot.dat')
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    full_command = f'transform {command}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Transform(), full_command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


def do_rotation(traj=None, command="", frame_indices=None, top=None):
    """Internal function to perform rotation (alias for rotate)"""
    return rotate(traj=traj, command=command, frame_indices=frame_indices, top=top)


def do_scaling(traj=None, command="", frame_indices=None, top=None):
    """Internal function to perform scaling (alias for scale)"""
    return scale(traj=traj, command=command, frame_indices=frame_indices, top=top)


def do_autoimage(traj=None, command="", frame_indices=None, top=None):
    """Internal function to perform autoimage (alias for autoimage)"""
    return autoimage(traj=traj, command=command, frame_indices=frame_indices, top=top)


@register_pmap
@super_dispatch()
def rotate_dihedral(traj=None, mask="", angle=0.0, frame_indices=None, top=None):
    """rotate around a dihedral angle

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        four-atom mask defining the dihedral
    angle : float, default 0.0
        rotation angle in degrees
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Rotate chi1 dihedral by 60 degrees
    >>> pt.rotate_dihedral(traj, ':2@N :2@CA :2@CB :2@CG', 60.0)
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'dihedral {mask} {angle}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Dihedral(), command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def set_dihedral(traj=None, mask="", angle=0.0, frame_indices=None, top=None):
    """set a dihedral angle to specific value

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        four-atom mask defining the dihedral
    angle : float, default 0.0
        target angle in degrees
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Set phi angle to -60 degrees
    >>> pt.set_dihedral(traj, ':2@C :3@N :3@CA :3@C', -60.0)
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'dihedral {mask} {angle} set'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Dihedral(), command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def replicate_cell(traj=None, command="", frame_indices=None, top=None):
    """replicate unit cell

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        replication command (e.g., 'nx 2 ny 2 nz 1')
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory with replicated coordinates

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Replicate 2x2x1
    >>> new_traj = pt.replicate_cell(traj, 'nx 2 ny 2 nz 1')
    """
    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    full_command = f'replicate {command}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Replicate(), full_command, top, dslist=dslist)

    # Create new trajectory with replicated coordinates
    new_frames = []
    for frame in traj:
        new_frame = frame.copy()
        action_list.compute_for_frame(new_frame)
        new_frames.append(new_frame)

    return Trajectory(xyz=np.array([f.xyz for f in new_frames]), top=top)


@register_pmap
@super_dispatch()
def randomize_ions(traj=None, mask="", around_mask="", distance=5.0,
                   frame_indices=None, top=None):
    """randomize positions of ions around a mask

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        ion mask to randomize
    around_mask : str
        mask to randomize around
    distance : float, default 5.0
        minimum distance from around_mask
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Randomize sodium ions around protein
    >>> pt.randomize_ions(traj, ':Na+', ':1-13', distance=8.0)
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'randomizeions {mask} around {around_mask} {distance}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_RandomizeIons(), command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def fiximagedbonds(traj=None, frame_indices=None, top=None):
    """fix imaged bonds across periodic boundaries

    Parameters
    ----------
    traj : Trajectory-like
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory (if traj is mutable) or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> pt.fiximagedbonds(traj)
    """
    _assert_mutable(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = 'fiximagedbonds'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_FixImagedBonds(), command, top, dslist=dslist)

    for frame in traj:
        action_list.compute_for_frame(frame)


@register_pmap
@super_dispatch()
def xtalsymm(traj=None, command="", frame_indices=None, top=None):
    """apply crystal symmetry operations

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        symmetry command options
    frame_indices : array-like, optional
    top : Topology, optional

    Returns
    -------
    Trajectory with symmetry applied

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Apply P21 symmetry
    >>> new_traj = pt.xtalsymm(traj, 'P21')
    """
    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    full_command = f'xtalsymm {command}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_XtalSymm(), full_command, top, dslist=dslist)

    # Create new trajectory with symmetry applied
    new_frames = []
    for frame in traj:
        new_frame = frame.copy()
        action_list.compute_for_frame(new_frame)
        new_frames.append(new_frame)

    return Trajectory(xyz=np.array([f.xyz for f in new_frames]), top=top)


# Export all functions
__all__ = [
    'align', 'superpose', 'center', 'translate', 'rotate', 'scale',
    'autoimage', 'image', 'transform', 'do_scaling', 'do_rotation',
    'do_autoimage', 'rotate_dihedral', 'set_dihedral', 'replicate_cell',
    'randomize_ions', 'fiximagedbonds', 'xtalsymm'
]