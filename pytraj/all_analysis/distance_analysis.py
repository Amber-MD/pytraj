"""
Distance-based analysis functions
"""
import numpy as np
from functools import partial
from typing import Union

from ..utils.get_common_objects import (
    get_topology, get_fiterator, get_data_from_dtype,
    get_list_of_commands, get_reference
)
from ..utils import ensure_not_none_or_string
from ..utils.decorators import register_pmap
from ..utils.convert import array_to_cpptraj_atommask
from ..datasets.datasetlist import DatasetList
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from ..trajectory.shared_methods import iterframe_master
from ..analysis.c_action import c_action
from ..analysis.c_action.actionlist import ActionList
from .base_classes import DatasetType


def pair_distance(p1, p2):
    """Calculate distance between two 3D points"""
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)


def in_voxel(voxel_cntr, xyz, delta):
    """Check if xyz coordinate is within voxel boundaries"""
    return (xyz[0] >= voxel_cntr[0] - delta and xyz[0] <= voxel_cntr[0] +
            delta) and (xyz[1] >= voxel_cntr[1] - delta
                        and xyz[1] <= voxel_cntr[1] + delta) and (
                            xyz[2] >= voxel_cntr[2] - delta
                            and xyz[2] <= voxel_cntr[2] + delta)


def _calculate_distance(traj, int_2darr: np.ndarray, n_frames: int, dtype: str) -> Union[np.ndarray, DatasetList]:
    """Internal function to calculate distances using integer arrays"""
    if int_2darr.ndim == 1:
        int_2darr = np.atleast_2d(int_2darr)

    arr = np.empty([n_frames, len(int_2darr)])

    for idx, frame in enumerate(iterframe_master(traj)):
        arr[idx] = frame._distance(int_2darr)

    arr = arr.T
    if dtype == 'ndarray':
        return arr
    else:
        dslist = DatasetList({'distance': arr})
        return get_data_from_dtype(dslist, dtype)


@register_pmap
def distance(traj=None,
             mask="",
             frame_indices=None,
             dtype='ndarray',
             top=None,
             image=False,
             n_frames=None):
    """compute distance between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or a list of string or a 2D array-like of integers.
        If `mask` is a 2D-array, the `image` option is always `False`.
        In this case, make sure to `autoimage` your trajectory before
        calling `distance`.
    frame_indices : array-like, optional, default None
    dtype : return type, default 'ndarray'
    top : Topology, optional
    image : bool, default False
    n_frames : int, optional, default None
        only need to provide n_frames if ``traj`` does not have this info

    Returns
    -------
    1D ndarray if mask is a string
    2D ndarray, shape (n_atom_pairs, n_frames) if mask is a list of strings or an array

    Notes
    -----
    Be careful with Topology. If your topology has Box info but your traj does not, you
    would get weird output ([0.0, ...]). Make sure to use `image=False` in this method or
    set_nobox for Topology.

    Examples
    --------
    >>> import pytraj as pt
    >>> # calculate distance for two atoms, using amber mask
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> dist = pt.distance(traj, '@1 @3')

    >>> # calculate distance for two groups of atoms, using amber mask
    >>> dist = pt.distance(traj, '@1,37,8 @2,4,6')

    >>> # calculate distance between two residues, using amber mask
    >>> dist = pt.distance(traj, ':1 :10')

    >>> # calculate multiple distances between two residues, using amber mask
    >>> # distance between residue 1 and 10, distance between residue 3 and 20
    >>> # (when using atom string mask, index starts from 1)
    >>> dist = pt.distance(traj, [':1 :10', ':3 :20'])

    >>> # calculate distance for a series of atoms, using array for atom mask
    >>> # distance between atom 1 and 5, distance between atom 4 and 10 (index starts from 0)
    >>> dist = pt.distance(traj, [[1, 5], [4, 10]])
    """
    ensure_not_none_or_string(traj)
    command = mask

    traj = get_fiterator(traj, frame_indices)
    topology = get_topology(traj, top)
    noimage_str = 'noimage' if not image or (hasattr(traj, 'crdinfo') and not traj.crdinfo['has_box'] and topology.has_box()) else ''

    command_array = np.asarray(command)

    if isinstance(command, (list, tuple, str, np.ndarray, int)):
        if 'int' in command_array.dtype.name:
            integer_array = command_array
            frame_count = traj.n_frames if n_frames is None else n_frames
            return _calculate_distance(traj, integer_array, frame_count, dtype)
        else:
            command_list = get_list_of_commands(command) if not isinstance(command, np.ndarray) else command
            cpptraj_action_datasets = CpptrajDatasetList()
            action_list = ActionList()

            for command in command_list:
                if noimage_str:
                    command = ' '.join((command, noimage_str))
                action_list.add(c_action.Action_Distance(), command, topology, dslist=cpptraj_action_datasets)

            action_list.compute(traj)
            return get_data_from_dtype(cpptraj_action_datasets, dtype)
    else:
        raise ValueError(
            "command must be a string, a list/tuple of strings, or "
            "a numpy 2D array")


@register_pmap
def _distance_to_ref_or_point(traj=None,
                              mask="",
                              point=None,
                              ref=None,
                              dtype='ndarray'):
    """Internal function for distance to reference or point calculations"""
    if point and ref:
        raise ValueError("Must be either point or ref, not both")

    top = traj.top
    command = mask
    if not isinstance(command, np.ndarray):
        list_of_commands = get_list_of_commands(command)
    else:
        list_of_commands = command

    action_datasets = CpptrajDatasetList()

    if ref is not None:
        refname = 'myref'
        ref_top = ref.top or traj.top
        action_datasets.add(DatasetType.REFERENCE, name=refname)
        action_datasets[0].top = ref_top
        action_datasets[0].add_frame(ref)
    actlist = ActionList()

    for cm in list_of_commands:
        if point is not None:
            point_str = 'point ' + ' '.join(str(_) for _ in point)
            cm = ' '.join((cm, point_str))
        elif ref is not None:
            cm = ' '.join((cm, 'reference'))
        actlist.add(c_action.Action_Distance(), cm, top, dslist=action_datasets)

    for frame in traj:
        actlist.compute(frame)

    if ref is not None:
        # remove ref
        action_datasets._pop(0)
    return get_data_from_dtype(action_datasets, dtype)


distance_to_point = partial(_distance_to_ref_or_point, ref=None)
distance_to_point.__doc__ = """
Examples
--------
    pytraj.distance_to_point(traj, ':1', point=[0., 0., 0.]) # doctest: +SKIP
"""

distance_to_reference = partial(_distance_to_ref_or_point, point=None)
distance_to_reference.__doc__ = """
Examples
--------
    pytraj.distance_to_reference(traj, '@1 @1', ref=ref) # doctest: +SKIP
"""


def pairwise_distance(traj=None,
                      mask_1='',
                      mask_2='',
                      top=None,
                      dtype='ndarray',
                      frame_indices=None):
    '''compute pairwise distance between atoms in mask_1 and atoms in mask_2

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask_1: string or 1D array-like
    mask_2: string or 1D array-like
    ...

    Returns
    -------
    out_1 : numpy array, shape (n_frames, n_atom_1, n_atom_2)
    out_2 : atom pairs, shape=(n_atom_1, n_atom_2, 2)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> mat = pt.pairwise_distance(traj, '@CA', '@CB')

    Notes
    -----
    This method is only fast for small number of atoms.
    '''
    from itertools import product

    top_ = get_topology(traj, top)
    indices_1 = top_.select(mask_1) if isinstance(mask_1, str) else mask_1
    indices_2 = top_.select(mask_2) if isinstance(mask_2, str) else mask_2
    arr = np.array(list(product(indices_1, indices_2)))
    mat = distance(
        traj, mask=arr, dtype=dtype, top=top_, frame_indices=frame_indices)
    mat = mat.T
    return (mat.reshape(mat.shape[0], len(indices_1), len(indices_2)),
            arr.reshape(len(indices_1), len(indices_2), 2))


@register_pmap
def closest_atom(top=None, frame=None, point=(0, 0, 0), mask=""):
    """for a given xyz coordinate in a frame, find the closest atom

    Parameters
    ----------
    top: Topology object
    frame: Frame object
    point: Tuple containing 3 elements, for x y and z coordinate of point.
    Default point is the origin
    mask: string containing atoms to find, in atom mask syntax. Default is empty string
    which contains all the atoms

    Returns
    -------
    Index of atom closest to point in xyz coordinate space.

    Notes
    -----
    Topology needs to contain atoms that match the atom mask passed in, and frame needs to have
    xyz coordinates for all atoms. Point should be a tuple of length 3 with format (x, y, z)

    Examples
    --------
    >>> import pytraj as pt
    >>> # find closest atom to origin in a given trajectory
    >>> traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
    >>> frame = traj[0]
    >>> pt.closest_atom(traj.top, traj[0], (0, 0, 0))
    205
    """
    if (len(top.atom_indices(mask)) == 0):
        raise ValueError(
            "Please pass in a topology file with atoms that match the mask in it"
        )

    closest_dist = None
    closest_idx = None
    atoms = top.atom_indices(mask)
    for atm in atoms:
        coord = frame.atom(atm)
        if ((closest_dist is None)
                or (pair_distance(coord, point) < closest_dist)):
            closest_dist = pair_distance(coord, point)
            closest_idx = atm

    return closest_idx


@register_pmap
def count_in_voxel(traj=None, mask="", voxel_cntr=(0, 0, 0), voxel_size=5):
    """For a voxel with center xyz and size voxel_size, find atoms that match a given mask
    that are contained in that voxel over the course of a trajectory.

    This analysis command is meant to be used with GIST analysis to estimate water residence time
    for a given region. When running GIST on a trajectory, we get a list of voxels and their
    associated rotational/translational entropy. To analyze how long solvent molecules are
    residing in various regions of the macromolecular surface, we can compute residence times
    using the survival time correlation function. This can be done by plotting the number of
    solvent molecules that reside in the voxel at frame 1 and then plotting the decay in the
    number of those original solvent molecules as they diffuse out of the region and are
    replaced. This can provide insight into the behavior of solvent atoms closely interfacing
    with the macromolecule, since the dynamic behavior of these solvent atoms differs greatly
    from the bulk. Read more about residence time in 10.1021/jp020100m section 3.5

    Parameters
    ---------
    traj: Trajectory object with a loaded topology object
    mask: Mask with atoms that exist in the topology file
    voxel_cntr: xyz coordinates to define the center of the voxel
    voxel_size: height/length/width of voxel

    Returns
    -------
    List of lists, idx i contains list of atoms in the voxel at frame i

    Examples
    --------
    >>> pop = pt.count_in_voxel(tz2_traj[0:5], tz2_traj.top, "", xyz, 3)
    >>> print([len(i) for i in pop])

    >>> # calculate residence time for water molecules inside a GIST voxel of interest
    >>> tz2_traj = pt.datafiles.load_tz2_ortho()
    >>> wat_atoms = tz2_traj.top.select(":WAT")
    >>> gist_voxel = (35.26, 38.23, 1.66)
    >>> pop = pt.count_in_voxel(tz2_traj, tz2_traj.top, ":WAT@O", gist_voxel, 10)
    >>> #
    >>> # print water molecules in voxel at frame 0.
    >>> # For survival time correlation function,
    >>> # plot the number of these particular beginning frame water atoms that remain in the voxel
    >>> # throughout the course of the simulation (if the voxel is part of the bulk solvent,
    >>> # the function will decay very quickly, but if it is near the interface, then decay
    >>> # will be slower as waters will stay in enthalpically favorable points).
    >>> orig_frame_waters = set(pop[0])
    >>> #
    >>> # NOTE: Reader will need to modify this to also exclude waters that leave in an intermediate frame
    >>> # but then return later. (e.g. water 250 is in frame 1, but not in frame 2, then back in frame 3
    >>> # it does not count as a "surviving" water that has remained in the voxel.
    >>> survival_decay = [len(orig_frame_waters.intersection(set(pop[i]))) for i in range(len(pop))]
    >>> print(survival_decay)
    """
    lives_in_voxel = []
    population = traj.top.atom_indices(mask)
    delta = voxel_size / 2

    for frame in traj:
        frame_voxAtoms = []
        for atm in population:
            coord = frame.atom(atm)
            if (in_voxel(voxel_cntr, coord, delta)):
                frame_voxAtoms.append(atm)
        lives_in_voxel.append(frame_voxAtoms)

    return lives_in_voxel


@register_pmap
def mindist(traj=None,
            command="",
            top=None,
            dtype='ndarray',
            frame_indices=None):
    '''compute mindist

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.mindist(traj, '@CA @H')
    '''
    from ..analysis.c_action import do_action
    from ..utils.convert import array2d_to_cpptraj_maskgroup

    if not isinstance(command, str):
        command = array2d_to_cpptraj_maskgroup(command)
    command = "mindist " + command

    action_datasets, _ = do_action(traj, command, c_action.Action_NativeContacts)
    return get_data_from_dtype(action_datasets, dtype)[-1]