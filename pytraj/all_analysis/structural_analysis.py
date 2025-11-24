"""
Structural analysis functions (angles, dihedrals, RMSD, etc.)
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
    _calculate_angles_for_int_array, DatasetType
)


def _calculate_dihedrals_for_int_array(traj, integer_array, n_frames, dtype):
    """Internal function to calculate dihedrals using integer arrays"""
    if integer_array.ndim == 1:
        integer_array = np.atleast_2d(integer_array)

    arr = np.empty([n_frames, len(integer_array)])

    for idx, frame in enumerate(iterframe_master(traj)):
        arr[idx] = frame._dihedral(integer_array)

    arr = arr.T
    if dtype == 'ndarray':
        return arr
    else:
        dslist = DatasetList({'dihedral': arr})
        return get_data_from_dtype(dslist, dtype)


def _dihedral_res(traj, mask=(), resid=0, dtype='ndarray', top=None):
    """Calculate dihedral angles for a specific residue"""
    top_ = get_topology(traj, top)
    resname = top_.residue(resid).name

    dihedral_dict = {
        'phi': [f':{resid-1}@C', f':{resid}@N', f':{resid}@CA', f':{resid}@C'],
        'psi': [f':{resid}@N', f':{resid}@CA', f':{resid}@C', f':{resid+1}@N'],
        'chi1': None,
        'chi2': None,
        'chi3': None,
        'chi4': None,
        'chi5': None,
    }

    # Define chi angles based on residue type
    if resname in ['ARG', 'LYS', 'MET', 'GLU', 'GLN', 'PRO']:
        dihedral_dict['chi1'] = [f':{resid}@N', f':{resid}@CA', f':{resid}@CB', f':{resid}@CG']
    elif resname in ['ILE', 'THR', 'VAL']:
        dihedral_dict['chi1'] = [f':{resid}@N', f':{resid}@CA', f':{resid}@CB', f':{resid}@CG1']
    elif resname in ['SER']:
        dihedral_dict['chi1'] = [f':{resid}@N', f':{resid}@CA', f':{resid}@CB', f':{resid}@OG']
    elif resname in ['CYS']:
        dihedral_dict['chi1'] = [f':{resid}@N', f':{resid}@CA', f':{resid}@CB', f':{resid}@SG']

    # Add chi2, chi3, etc. for longer side chains
    if resname in ['ARG', 'LYS', 'GLU', 'GLN', 'MET']:
        dihedral_dict['chi2'] = [f':{resid}@CA', f':{resid}@CB', f':{resid}@CG', f':{resid}@CD']

    results = {}
    for angle_name in mask if mask else dihedral_dict.keys():
        if dihedral_dict.get(angle_name):
            mask_str = ' '.join(dihedral_dict[angle_name])
            results[angle_name] = dihedral(traj, mask=mask_str, dtype=dtype, top=top_)

    return results


@register_pmap
def angle(traj=None, mask="", frame_indices=None, dtype='ndarray',
          top=None, *args, **kwargs):
    """compute angle between two maskes

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask : str or a list of string or a 2D array-like of integers.
    frame_indices : array-like, optional, default None
    dtype : return type, default 'ndarray'
    top : Topology, optional

    Returns
    -------
    1D ndarray if mask is a string
    2D ndarray, shape (n_atom_triplets, n_frames) if mask is a list of strings or an array

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> angle_data = pt.angle(traj, '@1 @2 @3')
    >>> angle_data = pt.angle(traj, [':2@CA :2@CB :2@CG', ':3@CA :3@CB :3@CG'])
    """
    ensure_not_none_or_string(traj)
    command = mask

    traj = get_fiterator(traj, frame_indices)
    topology = get_topology(traj, top)

    command_array = np.asarray(command)

    if isinstance(command, (list, tuple, str, np.ndarray, int)):
        if 'int' in command_array.dtype.name:
            integer_array = command_array
            frame_count = traj.n_frames
            return _calculate_angles_for_int_array(traj, integer_array, frame_count, dtype)
        else:
            command_list = get_list_of_commands(command) if not isinstance(command, np.ndarray) else command
            cpptraj_action_datasets = CpptrajDatasetList()
            action_list = ActionList()

            for cmd in command_list:
                action_list.add(c_action.Action_Angle(), cmd, topology, dslist=cpptraj_action_datasets)

            action_list.compute(traj)
            return get_data_from_dtype(cpptraj_action_datasets, dtype)
    else:
        raise ValueError(
            "command must be a string, a list/tuple of strings, or "
            "a numpy 2D array")


@register_pmap
def dihedral(traj=None, mask="", top=None, dtype='ndarray',
            frame_indices=None, *args, **kwargs):
    """compute dihedral angle between two maskes

    Parameters
    ----------
    traj : Trajectory-like
    mask : str or a list of string or a 2D array-like of integers.
    top : Topology, optional, default=None
    dtype : return type, default 'ndarray'
    frame_indices : array-like, optional, default None

    Returns
    -------
    1D ndarray if mask is a string
    2D ndarray, shape (n_atom_quartets, n_frames) if mask is a list of strings or an array

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> phi = pt.dihedral(traj, ':2@C :3@N :3@CA :3@C')
    >>> psi = pt.dihedral(traj, ':2@N :2@CA :2@C :3@N')
    >>> phi_psi = pt.dihedral(traj, [':2@C :3@N :3@CA :3@C', ':2@N :2@CA :2@C :3@N'])
    """
    ensure_not_none_or_string(traj)
    command = mask

    traj = get_fiterator(traj, frame_indices)
    topology = get_topology(traj, top)

    command_array = np.asarray(command)

    if isinstance(command, (list, tuple, str, np.ndarray, int)):
        if 'int' in command_array.dtype.name:
            integer_array = command_array
            frame_count = traj.n_frames
            return _calculate_dihedrals_for_int_array(traj, integer_array, frame_count, dtype)
        else:
            command_list = get_list_of_commands(command) if not isinstance(command, np.ndarray) else command
            cpptraj_action_datasets = CpptrajDatasetList()
            action_list = ActionList()

            for cmd in command_list:
                action_list.add(c_action.Action_Dihedral(), cmd, topology, dslist=cpptraj_action_datasets)

            action_list.compute(traj)
            return get_data_from_dtype(cpptraj_action_datasets, dtype)
    else:
        raise ValueError(
            "command must be a string, a list/tuple of strings, or "
            "a numpy 2D array")


@register_pmap
@super_dispatch()
def rmsd(traj=None, mask='', ref=None, ref_mask='', frame_indices=None,
         mass=False, top=None, dtype='ndarray', nofit=False):
    """compute rmsd

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default ''
        Atom mask
    ref : {Frame, int}, default None (first frame)
        Reference frame, either a Frame or a frame index
    ref_mask : str, optional, default ''
        if not given, use `mask`
    frame_indices : {None, array-like}, default None
        if not None, compute RMSD for given frames
    mass : bool, default False
        if True, use mass-weighted
    top : Topology, optional, default None
    dtype : str, default 'ndarray'
        return data type
    nofit : bool, default False
        if True, don't perform fitting

    Returns
    -------
    1D ndarray, shape=(n_frames,)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.rmsd(traj, mask='@CA')
    >>> data = pt.rmsd(traj, mask='@CA', ref=0)
    >>> data = pt.rmsd(traj, mask='@CA', ref=traj[-3])
    """
    ensure_not_none_or_string(traj)

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
    command = f'rmsd {mask} reference'
    if ref_mask != mask:
        command += f' refmask {ref_mask}'
    if mass:
        command += ' mass'
    if nofit:
        command += ' nofit'

    # Setup datasets and reference
    dslist = CpptrajDatasetList()
    dslist.add(DatasetType.REFERENCE, name='myref')
    dslist[0].top = ref_frame.top or top
    dslist[0].add_frame(ref_frame)

    # Create action
    action_list = ActionList()
    action_list.add(c_action.Action_Rmsd(), command, top, dslist=dslist)
    action_list.compute(traj)

    # Remove reference from dataset list
    dslist._pop(0)
    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def rmsd_nofit(traj=None, mask='', ref=None, ref_mask='', frame_indices=None,
               mass=False, top=None, dtype='ndarray'):
    """compute rmsd without fitting (translation and rotation)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default ''
    ref : {Frame, int}, default None (first frame)
    ref_mask : str, optional, default ''
    frame_indices : {None, array-like}, default None
    mass : bool, default False
    top : Topology, optional, default None
    dtype : str, default 'ndarray'

    Returns
    -------
    1D ndarray, shape=(n_frames,)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.rmsd_nofit(traj, mask='@CA')
    """
    return rmsd(traj=traj, mask=mask, ref=ref, ref_mask=ref_mask,
                frame_indices=frame_indices, mass=mass, top=top,
                dtype=dtype, nofit=True)


@register_pmap
@super_dispatch()
def rmsd_perres(traj=None, mask='', ref=None, ref_mask='', frame_indices=None,
                mass=False, top=None, dtype='ndarray', resrange=None):
    """compute RMSD per residue

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default ''
    ref : {Frame, int}, default None (first frame)
    ref_mask : str, optional, default ''
    frame_indices : {None, array-like}, default None
    mass : bool, default False
    top : Topology, optional, default None
    dtype : str, default 'ndarray'
    resrange : str, optional
        Residue range, e.g. '1-10'

    Returns
    -------
    2D ndarray, shape=(n_residues, n_frames)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.rmsd_perres(traj, mask='@CA')
    """
    ensure_not_none_or_string(traj)

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
    command = f'rmsd {mask} reference perres'
    if ref_mask != mask:
        command += f' refmask {ref_mask}'
    if mass:
        command += ' mass'
    if resrange:
        command += f' range {resrange}'

    # Setup datasets and reference
    dslist = CpptrajDatasetList()
    dslist.add(DatasetType.REFERENCE, name='myref')
    dslist[0].top = ref_frame.top or top
    dslist[0].add_frame(ref_frame)

    # Create action
    action_list = ActionList()
    action_list.add(c_action.Action_Rmsd(), command, top, dslist=dslist)
    action_list.compute(traj)

    # Remove reference from dataset list
    dslist._pop(0)
    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def rmsf(traj=None, mask="", top=None, dtype='ndarray',
         frame_indices=None, options=''):
    """compute atomicfluct (RMSF)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    top : Topology, optional
    dtype : str, default 'ndarray'
        return data type
    frame_indices : array-like, optional, default None
        frame indices
    options : str, optional
        extra cpptraj options

    Returns
    -------
    1D ndarray

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.rmsf(traj, '@CA')
    """
    ensure_not_none_or_string(traj)
    command = f'atomicfluct {mask} {options}'.strip()

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_AtomicFluct(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


# Create alias
atomicfluct = rmsf


@register_pmap
@super_dispatch()
def native_contacts(traj=None, mask="", ref=None, distance=7.0,
                   frame_indices=None, top=None, dtype='ndarray',
                   options=''):
    """compute native contacts

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for native contact analysis
    ref : Frame or int, default None
        reference frame
    distance : float, default 7.0
        distance cutoff for native contacts
    frame_indices : array-like, optional
    top : Topology, optional
    dtype : str, default 'ndarray'
    options : str, optional
        additional cpptraj options

    Returns
    -------
    1D ndarray

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.native_contacts(traj, mask=':1-13@CA')
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    # Get reference
    if ref is None:
        ref_frame = traj[0]
    elif hasattr(ref, 'xyz'):
        ref_frame = ref
    else:
        ref_frame = traj[ref]

    # Build command
    command = f'nativecontacts {mask} distance {distance} reference {options}'.strip()

    # Setup datasets and reference
    dslist = CpptrajDatasetList()
    dslist.add(DatasetType.REFERENCE, name='myref')
    dslist[0].top = ref_frame.top or top
    dslist[0].add_frame(ref_frame)

    # Create action
    action_list = ActionList()
    action_list.add(c_action.Action_NativeContacts(), command, top, dslist=dslist)
    action_list.compute(traj)

    # Remove reference from dataset list
    dslist._pop(0)
    return get_data_from_dtype(dslist, dtype)


def bfactors(traj=None, mask='@CA', dtype='ndarray', top=None):
    """compute bfactor

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '@CA'
    dtype : str, default 'ndarray'
    top : Topology, optional

    Returns
    -------
    1D ndarray

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.bfactors(traj)
    """
    return rmsf(traj, mask=mask, dtype=dtype, top=top) * 8.0 * np.pi**2 / 3.0


def check_structure(traj=None, mask='', top=None, frame_indices=None,
                   dtype='dataset'):
    """check structure for unusual bond lengths, angles, etc.

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    top : Topology, optional
    frame_indices : array-like, optional
    dtype : str, default 'dataset'

    Returns
    -------
    DatasetList

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.check_structure(traj)
    """
    command = f'checkstructure {mask}'.strip()
    dslist, _ = do_action(traj, command, c_action.Action_CheckStructure,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def check_chirality(traj=None, mask='', top=None, frame_indices=None,
                   dtype='dataset'):
    """check chirality of residues

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    top : Topology, optional
    frame_indices : array-like, optional
    dtype : str, default 'dataset'

    Returns
    -------
    DatasetList

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.check_chirality(traj)
    """
    command = f'checkchirality {mask}'.strip()
    dslist, _ = do_action(traj, command, c_action.Action_CheckChirality,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Export all functions
__all__ = [
    'angle', 'dihedral', 'rmsd', 'rmsd_nofit', 'rmsd_perres',
    'native_contacts', 'rmsf', 'atomicfluct', 'bfactors',
    'check_structure', 'check_chirality', '_dihedral_res',
    '_calculate_dihedrals_for_int_array'
]