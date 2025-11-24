"""
Conformational analysis functions
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
def pucker(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
           method='cremer', amplitude=False, theta=False, phi=False):
    """calculate ring pucker parameters

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        ring atom mask (should be 5 or 6 atoms)
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    method : str, default 'cremer'
        pucker calculation method ('cremer' or 'altona')
    amplitude : bool, default False
        calculate pucker amplitude
    theta : bool, default False
        calculate theta angle
    phi : bool, default False
        calculate phi angle

    Returns
    -------
    DatasetList
        pucker parameters

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Analyze sugar ring pucker (5-membered ring)
    >>> pucker_data = pt.pucker(traj, ':1@C1\',C2\',C3\',C4\',O4\'')
    >>> # Calculate all pucker parameters
    >>> pucker_data = pt.pucker(traj, ':1@C1\',C2\',C3\',C4\',O4\'',
    ...                        amplitude=True, theta=True, phi=True)
    >>> # Use Altona method
    >>> pucker_data = pt.pucker(traj, ':1@C1\',C2\',C3\',C4\',O4\'', method='altona')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['pucker', mask, method]
    if amplitude:
        command_parts.append('amplitude')
    if theta:
        command_parts.append('theta')
    if phi:
        command_parts.append('phi')

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Pucker(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def multidihedral(traj=None, dihedral_list=None, frame_indices=None, top=None,
                 dtype='ndarray', range360=False):
    """calculate multiple dihedral angles

    Parameters
    ----------
    traj : Trajectory-like
    dihedral_list : list of str
        list of four-atom masks defining dihedrals
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    range360 : bool, default False
        if True, return angles in 0-360 range instead of -180 to 180

    Returns
    -------
    2D ndarray
        dihedral angles, shape (n_frames, n_dihedrals)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Define multiple backbone dihedrals
    >>> dihedrals = [':2@C :3@N :3@CA :3@C',  # phi
    ...              ':3@N :3@CA :3@C :4@N',  # psi
    ...              ':4@C :5@N :5@CA :5@C']  # phi
    >>> angles = pt.multidihedral(traj, dihedrals)
    >>> # Use 0-360 degree range
    >>> angles = pt.multidihedral(traj, dihedrals, range360=True)
    """
    ensure_not_none_or_string(traj)

    if dihedral_list is None:
        raise ValueError("dihedral_list is required")

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    # Calculate each dihedral separately and combine
    all_dihedrals = []

    for dihedral_mask in dihedral_list:
        command = f'dihedral {dihedral_mask}'
        if range360:
            command += ' range360'

        dslist = CpptrajDatasetList()
        action_list = ActionList()
        action_list.add(c_action.Action_Dihedral(), command, top, dslist=dslist)
        action_list.compute(traj)

        dihedral_data = get_data_from_dtype(dslist, 'ndarray')
        all_dihedrals.append(dihedral_data)

    # Combine all dihedrals
    result = np.column_stack(all_dihedrals)

    if dtype == 'ndarray':
        return result
    else:
        dslist = DatasetList()
        dslist.add('multidihedral', result)
        return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def permute_dihedrals(traj=None, dihedral_list=None, frame_indices=None,
                     top=None, check_for_bonds=True, interval=1):
    """permute dihedral angles to sample conformational space

    Parameters
    ----------
    traj : Trajectory-like
    dihedral_list : list of str
        list of four-atom masks defining dihedrals to permute
    frame_indices : array-like, optional, default None
    top : Topology, optional
    check_for_bonds : bool, default True
        check for bonds when permuting
    interval : int, default 1
        interval for permutation

    Returns
    -------
    Trajectory
        trajectory with permuted dihedrals

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Permute side chain dihedrals
    >>> dihedrals = [':5@N :5@CA :5@CB :5@CG']
    >>> perm_traj = pt.permute_dihedrals(traj, dihedrals)
    >>> # Don't check for bonds (faster)
    >>> perm_traj = pt.permute_dihedrals(traj, dihedrals, check_for_bonds=False)
    """
    ensure_not_none_or_string(traj)

    if dihedral_list is None:
        raise ValueError("dihedral_list is required")

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    # Build command
    command_parts = ['permutedihedrals']
    for dihedral in dihedral_list:
        command_parts.append(dihedral)

    command_parts.extend(['interval', str(interval)])
    if check_for_bonds:
        command_parts.append('checkbonds')

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_PermuteDihedrals(), command, top, dslist=dslist)

    # Apply to each frame
    for frame in traj:
        action_list.compute_for_frame(frame)

    return traj


@register_pmap
@super_dispatch()
def analyze_modes(traj=None, mask='', modes_file='', frame_indices=None,
                 top=None, dtype='ndarray', beg=1, end=10):
    """analyze normal modes or principal components

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for mode analysis
    modes_file : str
        file containing modes/eigenvectors
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    beg : int, default 1
        first mode to analyze
    end : int, default 10
        last mode to analyze

    Returns
    -------
    DatasetList
        mode analysis results

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Analyze first 5 modes
    >>> mode_data = pt.analyze_modes(traj, '@CA', 'modes.dat', beg=1, end=5)
    >>> # Analyze all backbone atoms
    >>> mode_data = pt.analyze_modes(traj, '@N,CA,C', 'modes.dat')
    """
    ensure_not_none_or_string(traj)

    if not modes_file:
        raise ValueError("modes_file is required")

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command = f'analyze modes {modes_file} {mask} beg {beg} end {end}'

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_AnalyzeModes(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


def ramachandran(traj=None, resrange=None, frame_indices=None, top=None,
                dtype='ndarray'):
    """calculate Ramachandran angles (phi, psi)

    Parameters
    ----------
    traj : Trajectory-like
    resrange : str or list, optional
        residue range for Ramachandran analysis
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'

    Returns
    -------
    3D ndarray
        phi and psi angles, shape (n_frames, n_residues, 2)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate for all residues
    >>> rama = pt.ramachandran(traj)
    >>> # Calculate for specific residues
    >>> rama = pt.ramachandran(traj, resrange='2-10')
    >>> # Calculate for residue list
    >>> rama = pt.ramachandran(traj, resrange=[2, 5, 8, 10])
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    if resrange is None:
        # Use all protein residues
        resrange = '1-' + str(top.n_residues)
    elif isinstance(resrange, list):
        resrange = ','.join(map(str, resrange))

    # Calculate phi and psi for each residue
    phi_angles = []
    psi_angles = []

    # Parse residue range
    if '-' in resrange:
        start, end = map(int, resrange.split('-'))
        residues = list(range(start, end + 1))
    else:
        residues = [int(x) for x in resrange.split(',')]

    for res in residues:
        # Phi angle: C(i-1) - N(i) - CA(i) - C(i)
        if res > 1:
            phi_mask = f':{res-1}@C :{res}@N :{res}@CA :{res}@C'
            phi = dihedral(traj, phi_mask, frame_indices=None, top=top, dtype='ndarray')
            phi_angles.append(phi)
        else:
            phi_angles.append(np.full(len(traj), np.nan))

        # Psi angle: N(i) - CA(i) - C(i) - N(i+1)
        if res < top.n_residues:
            psi_mask = f':{res}@N :{res}@CA :{res}@C :{res+1}@N'
            psi = dihedral(traj, psi_mask, frame_indices=None, top=top, dtype='ndarray')
            psi_angles.append(psi)
        else:
            psi_angles.append(np.full(len(traj), np.nan))

    # Combine phi and psi
    phi_array = np.column_stack(phi_angles)
    psi_array = np.column_stack(psi_angles)
    rama_array = np.stack([phi_array, psi_array], axis=2)

    if dtype == 'ndarray':
        return rama_array
    else:
        dslist = DatasetList()
        dslist.add('ramachandran', rama_array)
        return get_data_from_dtype(dslist, dtype)


def secondary_structure_content(traj=None, frame_indices=None, top=None,
                               dtype='ndarray', assign_method='dssp'):
    """calculate secondary structure content over time

    Parameters
    ----------
    traj : Trajectory-like
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    assign_method : str, default 'dssp'
        method for secondary structure assignment

    Returns
    -------
    2D ndarray
        secondary structure content, shape (n_frames, n_ss_types)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate secondary structure content
    >>> ss_content = pt.secondary_structure_content(traj)
    """
    ensure_not_none_or_string(traj)

    command = f'secstruct {assign_method}'

    dslist, _ = do_action(traj, command, c_action.Action_SecStruct,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Import dihedral function from structural_analysis if needed
from .structural_analysis import dihedral

# Export all functions
__all__ = [
    'pucker', 'multidihedral', 'permute_dihedrals', 'analyze_modes',
    'ramachandran', 'secondary_structure_content'
]