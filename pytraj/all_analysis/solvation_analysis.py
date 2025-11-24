"""
Solvation and dynamics analysis functions
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
def gist(traj=None, solute_mask='', solvent_mask=':WAT', frame_indices=None,
         top=None, dtype='dataset', gridcntr=None, griddim=None,
         gridspacn=0.5, temp=300.0):
    """Grid Inhomogeneous Solvation Theory analysis

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask : str
        solute atom mask
    solvent_mask : str, default ':WAT'
        solvent molecule mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    gridcntr : tuple, optional
        grid center coordinates (x, y, z)
    griddim : tuple, optional
        grid dimensions (nx, ny, nz)
    gridspacn : float, default 0.5
        grid spacing in Angstroms
    temp : float, default 300.0
        temperature in Kelvin

    Returns
    -------
    DatasetList
        GIST analysis results including thermodynamic quantities

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # GIST analysis around protein
    >>> gist_data = pt.gist(traj, ':1-13', ':WAT')
    >>> # Custom grid parameters
    >>> gist_data = pt.gist(traj, ':1-13', ':WAT',
    ...                     gridcntr=(0, 0, 0), griddim=(40, 40, 40))
    >>> # Different temperature
    >>> gist_data = pt.gist(traj, ':1-13', ':WAT', temp=298.0)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['gist', f'solute {solute_mask}', f'solvent {solvent_mask}',
                    f'gridspacn {gridspacn}', f'temp {temp}']

    if gridcntr is not None:
        command_parts.extend(['gridcntr', f'{gridcntr[0]}', f'{gridcntr[1]}', f'{gridcntr[2]}'])
    if griddim is not None:
        command_parts.extend(['griddim', f'{griddim[0]}', f'{griddim[1]}', f'{griddim[2]}'])

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_GIST(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def diffusion(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
              time=1.0, lower=1, upper=None, max_dist=None):
    """calculate diffusion coefficients

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for diffusion calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    time : float, default 1.0
        time step between frames in ps
    lower : int, default 1
        lower bound for MSD calculation
    upper : int, optional
        upper bound for MSD calculation
    max_dist : float, optional
        maximum distance for unwrapping

    Returns
    -------
    DatasetList
        diffusion coefficients and mean square displacement

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Diffusion of water molecules
    >>> diff_data = pt.diffusion(traj, ':WAT@O', time=0.5)
    >>> # Diffusion with custom parameters
    >>> diff_data = pt.diffusion(traj, ':WAT@O', time=1.0, upper=100)
    >>> # Ion diffusion
    >>> diff_data = pt.diffusion(traj, ':Na+', time=1.0)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['diffusion', mask, f'time {time}', f'lower {lower}']
    if upper is not None:
        command_parts.extend(['upper', str(upper)])
    if max_dist is not None:
        command_parts.extend(['maxdist', str(max_dist)])

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Diffusion(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def rotdif(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
           time=1.0, order=2, nvec=1000):
    """calculate rotational diffusion

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for rotational diffusion
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    time : float, default 1.0
        time step between frames in ps
    order : int, default 2
        order of Legendre polynomial
    nvec : int, default 1000
        number of random vectors for averaging

    Returns
    -------
    DatasetList
        rotational diffusion data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Rotational diffusion of water
    >>> rotdif_data = pt.rotdif(traj, ':WAT', time=0.5)
    >>> # Different order Legendre polynomial
    >>> rotdif_data = pt.rotdif(traj, ':WAT', time=1.0, order=1)
    >>> # More vectors for better averaging
    >>> rotdif_data = pt.rotdif(traj, ':WAT', time=1.0, nvec=5000)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'rotdif {mask} time {time} order {order} nvec {nvec}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_RotDif(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def ti(traj=None, mask1='', mask2='', frame_indices=None, top=None,
       dtype='ndarray', curve='linear', npoints=11):
    """thermodynamic integration analysis

    Parameters
    ----------
    traj : Trajectory-like
    mask1 : str
        first state mask
    mask2 : str
        second state mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    curve : str, default 'linear'
        integration curve type ('linear', 'polynomial')
    npoints : int, default 11
        number of integration points

    Returns
    -------
    DatasetList
        thermodynamic integration results

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # TI between two states
    >>> ti_data = pt.ti(traj, ':LIG', ':LIG_mod')
    >>> # Polynomial curve
    >>> ti_data = pt.ti(traj, ':LIG', ':LIG_mod', curve='polynomial')
    >>> # More integration points
    >>> ti_data = pt.ti(traj, ':LIG', ':LIG_mod', npoints=21)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'ti {mask1} {mask2} curve {curve} npoints {npoints}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_TI(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


def solvation_shell(traj=None, solute_mask='', solvent_mask=':WAT',
                   frame_indices=None, top=None, dtype='ndarray',
                   first_shell=3.5, second_shell=5.0):
    """analyze solvation shell structure

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask : str
        solute atom mask
    solvent_mask : str, default ':WAT'
        solvent molecule mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    first_shell : float, default 3.5
        first solvation shell distance cutoff
    second_shell : float, default 5.0
        second solvation shell distance cutoff

    Returns
    -------
    DatasetList
        solvation shell analysis

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Analyze water solvation around protein
    >>> shell_data = pt.solvation_shell(traj, ':1-13', ':WAT@O')
    >>> # Custom shell distances
    >>> shell_data = pt.solvation_shell(traj, ':1-13', ':WAT@O',
    ...                                first_shell=3.0, second_shell=6.0)
    """
    ensure_not_none_or_string(traj)

    command = (f'solvent {solute_mask} {solvent_mask} '
               f'first {first_shell} second {second_shell}')

    dslist, _ = do_action(traj, command, c_action.Action_Solvent,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def residence_time(traj=None, solute_mask='', solvent_mask=':WAT',
                  frame_indices=None, top=None, dtype='ndarray',
                  distance=3.5, time_window=10):
    """calculate solvent residence times

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask : str
        solute atom mask
    solvent_mask : str, default ':WAT'
        solvent molecule mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    distance : float, default 3.5
        distance cutoff for residence
    time_window : int, default 10
        time window for continuous residence

    Returns
    -------
    DatasetList
        residence time analysis

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Water residence times around protein
    >>> res_time = pt.residence_time(traj, ':1-13', ':WAT@O')
    >>> # Different parameters
    >>> res_time = pt.residence_time(traj, ':1-13', ':WAT@O',
    ...                             distance=4.0, time_window=20)
    """
    ensure_not_none_or_string(traj)

    command = (f'lifetime {solute_mask} {solvent_mask} '
               f'distance {distance} window {time_window}')

    dslist, _ = do_action(traj, command, c_action.Action_Lifetime,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Export all functions
__all__ = [
    'gist', 'diffusion', 'rotdif', 'ti', 'solvation_shell', 'residence_time'
]