"""
Geometry analysis functions (center of mass, radius of gyration, etc.)
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
def center_of_mass(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray'):
    """compute center of mass

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    2D ndarray, shape=(n_frames, 3)
        center of mass coordinates for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.center_of_mass(traj, ':1-13')
    >>> data = pt.center_of_mass(traj, '@CA')
    """
    ensure_not_none_or_string(traj)
    command = f'center {mask} mass origin'

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Center(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def center_of_geometry(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray'):
    """compute center of geometry (geometric center)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    2D ndarray, shape=(n_frames, 3)
        center of geometry coordinates for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.center_of_geometry(traj, ':1-13')
    >>> data = pt.center_of_geometry(traj, '@CA')
    """
    ensure_not_none_or_string(traj)
    command = f'center {mask} origin'

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Center(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def radgyr(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
           mass=False, tensor=False):
    """compute radius of gyration

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    mass : bool, default False
        if True, use mass-weighted calculation
    tensor : bool, default False
        if True, return tensor components

    Returns
    -------
    1D ndarray if tensor=False
    2D ndarray if tensor=True, shape=(n_frames, 6) for tensor components

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.radgyr(traj, '@CA')
    >>> data = pt.radgyr(traj, ':1-13', mass=True)
    >>> tensor_data = pt.radgyr(traj, '@CA', tensor=True)
    """
    ensure_not_none_or_string(traj)

    # Build command
    command_parts = ['radgyr', mask]
    if mass:
        command_parts.append('mass')
    if tensor:
        command_parts.append('tensor')

    command = ' '.join(command_parts)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Radgyr(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def radgyr_tensor(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray', mass=False):
    """compute radius of gyration tensor components

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    mass : bool, default False
        if True, use mass-weighted calculation

    Returns
    -------
    2D ndarray, shape=(n_frames, 6)
        tensor components: [Rg, Rgxx, Rgyy, Rgzz, Rgxy, Rgxz, Rgyz]

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.radgyr_tensor(traj, '@CA')
    """
    return radgyr(traj=traj, mask=mask, frame_indices=frame_indices,
                  top=top, dtype=dtype, mass=mass, tensor=True)


@register_pmap
@super_dispatch()
def principal_axes(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
                   mass=False, dorotation=False):
    """compute principal axes of inertia

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    mass : bool, default False
        if True, use mass-weighted calculation
    dorotation : bool, default False
        if True, perform rotation to align with principal axes

    Returns
    -------
    DatasetList containing principal axes data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.principal_axes(traj, '@CA')
    >>> data = pt.principal_axes(traj, ':1-13', mass=True)
    """
    ensure_not_none_or_string(traj)

    # Build command
    command_parts = ['principal', mask]
    if mass:
        command_parts.append('mass')
    if dorotation:
        command_parts.append('dorotation')

    command = ' '.join(command_parts)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Principal(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def align_principal_axis(traj=None, mask='', frame_indices=None, top=None,
                        mass=False, axis='z'):
    """align trajectory along principal axis

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for alignment
    frame_indices : array-like, optional, default None
    top : Topology, optional
    mass : bool, default False
        if True, use mass-weighted calculation
    axis : str, default 'z'
        axis to align along ('x', 'y', or 'z')

    Returns
    -------
    Trajectory
        aligned trajectory

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> aligned_traj = pt.align_principal_axis(traj, '@CA')
    """
    ensure_not_none_or_string(traj)
    from ..trajectory import Trajectory

    # Build command
    command_parts = ['principal', mask]
    if mass:
        command_parts.append('mass')
    command_parts.extend(['dorotation', f'axis{axis}'])

    command = ' '.join(command_parts)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    # Create a copy of trajectory to modify
    aligned_traj = Trajectory(xyz=traj.xyz.copy(), top=top)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Principal(), command, top, dslist=dslist)

    # Apply transformation to each frame
    for frame in aligned_traj:
        action_list.compute_for_frame(frame)

    return aligned_traj


def moments_of_inertia(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray'):
    """compute moments of inertia

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    DatasetList containing moments of inertia

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.moments_of_inertia(traj, '@CA')
    """
    ensure_not_none_or_string(traj)
    command = f'inertialtensor {mask}'

    dslist, _ = do_action(traj, command, c_action.Action_InertialTensor,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def gyrate(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
           mass=False):
    """compute radius of gyration (alias for radgyr)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    mass : bool, default False
        if True, use mass-weighted calculation

    Returns
    -------
    1D ndarray

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.gyrate(traj, '@CA')
    """
    return radgyr(traj=traj, mask=mask, frame_indices=frame_indices,
                  top=top, dtype=dtype, mass=mass)


def bounds(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray'):
    """compute coordinate bounds (min/max for x, y, z)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    DatasetList containing bounds information

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.bounds(traj, '@CA')
    """
    ensure_not_none_or_string(traj)
    command = f'bounds {mask}'

    dslist, _ = do_action(traj, command, c_action.Action_Bounds,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def vector_analysis(traj=None, mask1='', mask2='', frame_indices=None,
                   top=None, dtype='ndarray', center=False,
                   principal_axes_reference=False):
    """analyze vectors between two masks

    Parameters
    ----------
    traj : Trajectory-like
    mask1 : str
        first atom mask (vector origin)
    mask2 : str
        second atom mask (vector end)
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    center : bool, default False
        if True, center vectors
    principal_axes_reference : bool, default False
        if True, align to principal axes

    Returns
    -------
    DatasetList containing vector data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.vector_analysis(traj, ':1@CA', ':2@CA')
    """
    ensure_not_none_or_string(traj)

    command_parts = ['vector', mask1, mask2]
    if center:
        command_parts.append('center')
    if principal_axes_reference:
        command_parts.append('principal')

    command = ' '.join(command_parts)

    dslist, _ = do_action(traj, command, c_action.Action_Vector,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Export all functions
__all__ = [
    'center_of_mass', 'center_of_geometry', 'radgyr', 'radgyr_tensor',
    'principal_axes', 'align_principal_axis', 'moments_of_inertia',
    'gyrate', 'bounds', 'vector_analysis'
]