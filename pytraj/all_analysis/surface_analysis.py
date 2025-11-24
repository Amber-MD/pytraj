"""
Surface and volume analysis functions
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
def surf(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
         probe=1.4):
    """compute solvent accessible surface area

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for surface calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    probe : float, default 1.4
        probe radius in Angstroms

    Returns
    -------
    1D ndarray
        surface area values for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate surface area for entire system
    >>> data = pt.surf(traj)
    >>> # Calculate surface area for protein only
    >>> data = pt.surf(traj, ':1-13')
    >>> # Use different probe radius
    >>> data = pt.surf(traj, ':1-13', probe=1.2)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'surf {mask}'
    if probe != 1.4:
        command += f' probe {probe}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Surf(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def molsurf(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
            probe=1.4):
    """compute molecular surface area

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for surface calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    probe : float, default 1.4
        probe radius in Angstroms

    Returns
    -------
    1D ndarray
        molecular surface area values for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate molecular surface area
    >>> data = pt.molsurf(traj, ':1-13')
    >>> # Use custom probe radius
    >>> data = pt.molsurf(traj, ':1-13', probe=1.2)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'molsurf {mask}'
    if probe != 1.4:
        command += f' probe {probe}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Molsurf(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def volume(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray'):
    """compute volume of atoms in mask

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for volume calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    1D ndarray
        volume values for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate volume of protein
    >>> data = pt.volume(traj, ':1-13')
    >>> # Calculate volume of specific residues
    >>> data = pt.volume(traj, ':1,5,10')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'volume {mask}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Volume(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def volmap(traj=None, mask='', frame_indices=None, top=None, dtype='dataset',
           dx=0.5, buffer=3.0, centermask=''):
    """create a volume map (3D grid) of atom density

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for volume map
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    dx : float, default 0.5
        grid spacing in Angstroms
    buffer : float, default 3.0
        buffer around atoms in Angstroms
    centermask : str, optional
        mask to center the grid on

    Returns
    -------
    DatasetList
        volume map data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Create volume map of water
    >>> data = pt.volmap(traj, ':WAT@O')
    >>> # Create high-resolution map
    >>> data = pt.volmap(traj, ':WAT@O', dx=0.25)
    >>> # Center on protein
    >>> data = pt.volmap(traj, ':WAT@O', centermask=':1-13')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['volmap', mask, f'dx {dx}', f'buffer {buffer}']
    if centermask:
        command_parts.extend(['centermask', centermask])

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Volmap(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def density(traj=None, mask='', frame_indices=None, top=None, dtype='dataset',
            delta=0.25, x=None, y=None, z=None):
    """calculate density distribution on a 3D grid

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for density calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    delta : float, default 0.25
        grid spacing
    x : tuple, optional
        x-axis range (min, max)
    y : tuple, optional
        y-axis range (min, max)
    z : tuple, optional
        z-axis range (min, max)

    Returns
    -------
    DatasetList
        density data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate water density
    >>> data = pt.density(traj, ':WAT@O')
    >>> # Calculate with custom grid spacing
    >>> data = pt.density(traj, ':WAT@O', delta=0.5)
    >>> # Calculate in specific region
    >>> data = pt.density(traj, ':WAT@O', x=(-10, 10), y=(-10, 10), z=(-10, 10))
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['density', mask, f'delta {delta}']

    if x is not None:
        command_parts.extend(['x', f'{x[0]}', f'{x[1]}'])
    if y is not None:
        command_parts.extend(['y', f'{y[0]}', f'{y[1]}'])
    if z is not None:
        command_parts.extend(['z', f'{z[0]}', f'{z[1]}'])

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Density(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def watershell(traj=None, solute_mask='', solvent_mask=':WAT',
               frame_indices=None, top=None, dtype='ndarray',
               lower=0.0, upper=3.4):
    """calculate number of solvent molecules in shell around solute

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask : str
        mask for solute atoms
    solvent_mask : str, default ':WAT'
        mask for solvent molecules
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    lower : float, default 0.0
        lower distance cutoff in Angstroms
    upper : float, default 3.4
        upper distance cutoff in Angstroms

    Returns
    -------
    1D ndarray
        number of solvent molecules in shell for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate water molecules around protein
    >>> data = pt.watershell(traj, ':1-13', ':WAT')
    >>> # Calculate with custom distance range
    >>> data = pt.watershell(traj, ':1-13', ':WAT', lower=2.0, upper=5.0)
    >>> # Calculate around specific residue
    >>> data = pt.watershell(traj, ':1@CA', ':WAT')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'watershell {solute_mask} {solvent_mask} lower {lower} upper {upper}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Watershell(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def grid(traj=None, mask='', frame_indices=None, top=None, dtype='dataset',
         spacing=1.0, box=None, center=None, pdb=False):
    """create a 3D grid for various analyses

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for grid calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    spacing : float, default 1.0
        grid spacing in Angstroms
    box : tuple, optional
        box dimensions (x, y, z)
    center : tuple, optional
        center coordinates (x, y, z)
    pdb : bool, default False
        output in PDB format

    Returns
    -------
    DatasetList
        grid data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Create basic grid
    >>> data = pt.grid(traj, ':WAT@O')
    >>> # Create with custom spacing
    >>> data = pt.grid(traj, ':WAT@O', spacing=0.5)
    >>> # Create with custom box
    >>> data = pt.grid(traj, ':WAT@O', box=(20, 20, 20))
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['grid', mask, f'spacing {spacing}']

    if box is not None:
        command_parts.extend(['box', f'{box[0]}', f'{box[1]}', f'{box[2]}'])
    if center is not None:
        command_parts.extend(['center', f'{center[0]}', f'{center[1]}', f'{center[2]}'])
    if pdb:
        command_parts.append('pdb')

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Grid(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


def _grid(traj=None, mask='', frame_indices=None, top=None, dtype='dataset',
          **kwargs):
    """Internal grid function with additional options

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for grid calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    **kwargs : dict
        additional grid options

    Returns
    -------
    DatasetList
        grid data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt._grid(traj, ':WAT@O', spacing=0.25, normalize=True)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['grid', mask]

    # Add all keyword arguments to command
    for key, value in kwargs.items():
        if isinstance(value, bool):
            if value:
                command_parts.append(key)
        else:
            command_parts.extend([key, str(value)])

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Grid(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


def cavity_volume(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
                  probe=1.4):
    """calculate cavity volume within a mask

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask defining the cavity boundary
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    probe : float, default 1.4
        probe radius for cavity detection

    Returns
    -------
    1D ndarray
        cavity volume for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate cavity volume in protein
    >>> data = pt.cavity_volume(traj, ':1-13')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'volume {mask} cavity probe {probe}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Volume(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


def surface_tension(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray'):
    """calculate surface tension

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for surface tension calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    1D ndarray
        surface tension values for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.surface_tension(traj, ':WAT')
    """
    ensure_not_none_or_string(traj)

    command = f'surftens {mask}'
    dslist, _ = do_action(traj, command, c_action.Action_SurfTens,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Export all functions
__all__ = [
    'surf', 'molsurf', 'volume', 'volmap', 'density', 'watershell',
    'grid', '_grid', 'cavity_volume', 'surface_tension'
]