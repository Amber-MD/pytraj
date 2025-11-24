"""
Energy analysis functions
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
def esander(traj=None, prmtop=None, igb=1, saltcon=0.0, cut=8.0,
           frame_indices=None, dtype='ndarray', ntb=0, dielc=1.0):
    """calculate sander energies using AMBER

    Parameters
    ----------
    traj : Trajectory-like
    prmtop : str or Topology
        parameter/topology file
    igb : int, default 1
        generalized Born model (0=vacuum, 1=HCT, 2=OBC1, 5=OBC2, 7=GBn, 8=GBn2)
    saltcon : float, default 0.0
        salt concentration for GB calculations
    cut : float, default 8.0
        nonbonded cutoff distance
    frame_indices : array-like, optional, default None
    dtype : str, default 'ndarray'
    ntb : int, default 0
        periodic boundary conditions (0=none, 1=constant volume, 2=constant pressure)
    dielc : float, default 1.0
        dielectric constant

    Returns
    -------
    DatasetList
        energy components (total, kinetic, potential, etc.)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate energies using GB model
    >>> energies = pt.esander(traj, traj.top, igb=1)
    >>> # Calculate vacuum energies
    >>> energies = pt.esander(traj, traj.top, igb=0)
    >>> # Use custom parameters
    >>> energies = pt.esander(traj, traj.top, igb=5, saltcon=0.15, cut=12.0)

    Notes
    -----
    Requires AMBER installation with sander. The trajectory must be compatible
    with the provided topology file.
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)

    # Handle topology
    if isinstance(prmtop, str):
        top = get_topology(traj, None)
        prmtop_file = prmtop
    else:
        top = get_topology(traj, prmtop)
        prmtop_file = prmtop.filename if hasattr(prmtop, 'filename') else None

    if prmtop_file is None:
        raise ValueError("prmtop file path is required for esander")

    # Build command
    command_parts = [
        'esander',
        f'prmtop {prmtop_file}',
        f'igb {igb}',
        f'saltcon {saltcon}',
        f'cut {cut}',
        f'ntb {ntb}',
        f'dielc {dielc}'
    ]

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Esander(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def lie(traj=None, mask1='', mask2='', frame_indices=None, top=None,
        dtype='ndarray', temp=300.0, prmtop=None):
    """calculate Linear Interaction Energy (LIE)

    Parameters
    ----------
    traj : Trajectory-like
    mask1 : str
        first mask (usually ligand)
    mask2 : str
        second mask (usually environment/protein+solvent)
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    temp : float, default 300.0
        temperature in Kelvin
    prmtop : str or Topology, optional
        parameter/topology file for energy calculations

    Returns
    -------
    DatasetList
        LIE energy components

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate LIE between ligand and environment
    >>> lie_data = pt.lie(traj, ':LIG', ':1-13,WAT', temp=298.0)
    >>> # Use with specific topology
    >>> lie_data = pt.lie(traj, ':LIG', ':1-13,WAT', prmtop='system.prmtop')

    Notes
    -----
    Linear Interaction Energy method for estimating binding free energies.
    Requires proper parameterization and sufficient sampling.
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    # Build command
    command_parts = ['lie', mask1, mask2, f'temp {temp}']

    if prmtop is not None:
        if isinstance(prmtop, str):
            command_parts.extend(['prmtop', prmtop])
        elif hasattr(prmtop, 'filename'):
            command_parts.extend(['prmtop', prmtop.filename])

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_LIE(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


def interaction_energy(traj=None, mask1='', mask2='', frame_indices=None,
                      top=None, dtype='ndarray', cut=8.0, elec=True,
                      vdw=True, nb_only=True):
    """calculate interaction energy between two masks

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
    cut : float, default 8.0
        cutoff distance for interactions
    elec : bool, default True
        include electrostatic interactions
    vdw : bool, default True
        include van der Waals interactions
    nb_only : bool, default True
        only calculate non-bonded interactions

    Returns
    -------
    1D ndarray
        interaction energies for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate protein-ligand interaction energy
    >>> ie = pt.interaction_energy(traj, ':1-13', ':LIG')
    >>> # Only electrostatic interactions
    >>> ie_elec = pt.interaction_energy(traj, ':1-13', ':LIG', vdw=False)
    >>> # Use longer cutoff
    >>> ie = pt.interaction_energy(traj, ':1-13', ':LIG', cut=12.0)
    """
    ensure_not_none_or_string(traj)

    # Build command
    command_parts = ['pairwise', mask1, mask2, f'cut {cut}']

    if not elec:
        command_parts.append('noelec')
    if not vdw:
        command_parts.append('novdw')
    if nb_only:
        command_parts.append('nb')

    command = ' '.join(command_parts)

    dslist, _ = do_action(traj, command, c_action.Action_Pairwise,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def decomp_interaction_energy(traj=None, mask1='', mask2='', frame_indices=None,
                             top=None, dtype='dataset', cut=8.0, byres=False):
    """decompose interaction energy by residue or atom

    Parameters
    ----------
    traj : Trajectory-like
    mask1 : str
        first atom mask
    mask2 : str
        second atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    cut : float, default 8.0
        cutoff distance for interactions
    byres : bool, default False
        decompose by residue instead of atom

    Returns
    -------
    DatasetList
        decomposed interaction energies

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Decompose by residue
    >>> decomp = pt.decomp_interaction_energy(traj, ':1-13', ':LIG', byres=True)
    >>> # Decompose by atom
    >>> decomp = pt.decomp_interaction_energy(traj, ':1-13', ':LIG')
    """
    ensure_not_none_or_string(traj)

    # Build command
    command_parts = ['pairwise', mask1, mask2, f'cut {cut}', 'decomp']

    if byres:
        command_parts.append('byres')

    command = ' '.join(command_parts)

    dslist, _ = do_action(traj, command, c_action.Action_Pairwise,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def energy_decomposition(traj=None, mask='', frame_indices=None, top=None,
                        dtype='dataset', prmtop=None, byres=False,
                        idecomp=1, dec_verbose=0):
    """perform energy decomposition analysis

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for decomposition
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    prmtop : str or Topology, optional
        parameter file
    byres : bool, default False
        decompose by residue
    idecomp : int, default 1
        decomposition method (1=per-residue, 2=per-residue pair, 3=pairwise, 4=all)
    dec_verbose : int, default 0
        verbosity level

    Returns
    -------
    DatasetList
        energy decomposition data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Per-residue decomposition
    >>> decomp = pt.energy_decomposition(traj, ':1-13', idecomp=1)
    >>> # Pairwise decomposition
    >>> decomp = pt.energy_decomposition(traj, ':1-13', idecomp=3)
    """
    ensure_not_none_or_string(traj)

    # Build command
    command_parts = ['energydecomp', mask, f'idecomp {idecomp}',
                    f'dec_verbose {dec_verbose}']

    if prmtop is not None:
        if isinstance(prmtop, str):
            command_parts.extend(['prmtop', prmtop])
        elif hasattr(prmtop, 'filename'):
            command_parts.extend(['prmtop', prmtop.filename])

    if byres:
        command_parts.append('byres')

    command = ' '.join(command_parts)

    dslist, _ = do_action(traj, command, c_action.Action_EnergyDecomp,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def potential_energy(traj=None, mask='', frame_indices=None, top=None,
                    dtype='ndarray', prmtop=None, cut=8.0):
    """calculate potential energy for atoms in mask

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for energy calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    prmtop : str or Topology, optional
        parameter file
    cut : float, default 8.0
        nonbonded cutoff

    Returns
    -------
    1D ndarray
        potential energies for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate potential energy for entire system
    >>> pe = pt.potential_energy(traj, ':*')
    >>> # Calculate for protein only
    >>> pe = pt.potential_energy(traj, ':1-13')
    """
    ensure_not_none_or_string(traj)

    # Build command
    command_parts = ['energy', mask, f'cut {cut}']

    if prmtop is not None:
        if isinstance(prmtop, str):
            command_parts.extend(['prmtop', prmtop])
        elif hasattr(prmtop, 'filename'):
            command_parts.extend(['prmtop', prmtop.filename])

    command = ' '.join(command_parts)

    dslist, _ = do_action(traj, command, c_action.Action_Energy,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def kinetic_energy(traj=None, mask='', frame_indices=None, top=None,
                  dtype='ndarray', temp=300.0):
    """calculate kinetic energy from velocities

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for kinetic energy calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    temp : float, default 300.0
        temperature for Maxwell-Boltzmann distribution

    Returns
    -------
    1D ndarray
        kinetic energies for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate kinetic energy (requires velocity data)
    >>> ke = pt.kinetic_energy(traj, '@CA')
    >>> # Use different temperature
    >>> ke = pt.kinetic_energy(traj, '@CA', temp=298.0)

    Notes
    -----
    Requires velocity information in the trajectory. If velocities are not
    present, they can be estimated from coordinate differences.
    """
    ensure_not_none_or_string(traj)

    command = f'kineticenergy {mask} temp {temp}'

    dslist, _ = do_action(traj, command, c_action.Action_KineticEnergy,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Export all functions
__all__ = [
    'esander', 'lie', 'interaction_energy', 'decomp_interaction_energy',
    'energy_decomposition', 'potential_energy', 'kinetic_energy'
]