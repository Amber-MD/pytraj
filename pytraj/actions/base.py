"""
Base utilities and common imports for pytraj actions
"""
from __future__ import absolute_import
import numpy as np
try:
    from enum import StrEnum
except ImportError:
    from enum import Enum
    class StrEnum(str, Enum):
        """
        Enum where members are also (and must be) strs.
        """
        def __new__(cls, value):
            member = str.__new__(cls, value)  # Create a new instance of str with the given value
            member._value_ = value  # Set the _value_ attribute to the given value
            return member

from typing import Any, Callable, List, Union
from functools import partial

from ..utils.get_common_objects import (
    get_topology,
    get_data_from_dtype,
    get_list_of_commands,
    get_reference,
    get_fiterator,
    super_dispatch, )
from ..utils import ensure_not_none_or_string
from ..utils import is_int
from ..utils.context import tempfolder
from ..utils.context import capture_stdout
from ..utils.convert import array_to_cpptraj_atommask
from ..utils.convert import array2d_to_cpptraj_maskgroup
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from ..datasets.datasetlist import DatasetList
from ..trajectory.shared_methods import iterframe_master
from ..trajectory.frame import Frame
from ..trajectory.trajectory import Trajectory
from ..trajectory.trajectory_iterator import TrajectoryIterator
from ..utils.decorators import register_pmap, register_openmp
from ..analysis.c_action import c_action
from ..analysis.c_action import do_action
from ..analysis.c_analysis import c_analysis
from ..analysis.c_action.actionlist import ActionList
from ..topology.topology import Topology
from ..builder.build import make_structure
from ..analysis.rmsd import (
    rotation_matrix,
    pairwise_rmsd,
    rmsd_perres,
    rmsd_nofit,
    rmsd,
    symmrmsd,
    distance_rmsd, )
from ..analysis.energy_analysis import (
    esander,
    lie, )
from ..analysis import (
    matrix,
    vector,
    nmr,
    dssp_analysis,
    hbond_analysis,
    energy_analysis, )

from ..core.c_core import CpptrajState, Command
from ..utils.command_builder import CommandBuilder


class DatasetType(StrEnum):
    COORDS = 'coords'
    REFERENCE = 'reference'
    REFERENCE_FRAME = 'ref_frame'
    TOPOLOGY = 'topology'
    VECTOR = 'vector'
    DOUBLE = 'double'
    MODES = 'modes'
    XYMESH = 'xymesh'
    MATRIX3x3 = 'matrix3x3'
    MATRIX_DBL = 'matrix_dbl'


def _assert_mutable(trajiter):
    """if pytraj.Trajectory : do nothing
    if pytraj.TrajReadOnly: copy to Trajectory
    """
    if isinstance(trajiter, TrajectoryIterator):
        return Trajectory(trajiter)
    else:
        return trajiter


def in_voxel(voxel_cntr, xyz, delta):
    x0, y0, z0 = voxel_cntr
    x, y, z = xyz

    return np.abs(x - x0) <= delta and np.abs(y - y0) <= delta and np.abs(z - z0) <= delta


def count_in_voxel(traj=None, mask="", voxel_cntr=(0, 0, 0), voxel_size=5):
    """count number of atoms in given voxel

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        mask
    voxel_cntr : tuple of floats, default (0, 0, 0)
        box center
    voxel_size : float, default 5.0
        box edge

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> count_in_voxel(traj, mask='@CA', voxel_cntr=(10, 20, 30), voxel_size=5)
    array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    >>> count_in_voxel(traj, mask='@CA', voxel_cntr=(10, 20, 30), voxel_size=50)
    array([12, 12, 12, 12, 12, 12, 12, 12, 12, 12])
    """
    half_size = voxel_size / 2.
    atom_indices = traj.top.select(mask)
    count_arr = np.empty(traj.n_frames, dtype='i')

    for n, frame in enumerate(traj):
        xyz_array = frame.xyz[atom_indices]
        count = 0
        for xyz in xyz_array:
            if in_voxel(voxel_cntr, xyz, half_size):
                count += 1
        count_arr[n] = count

    return count_arr


def pair_distance(p1, p2):
    '''
    Distance between two vectors
    '''
    return np.sqrt(np.sum((p1 - p2)**2))


def closest_atom(top=None, frame=None, point=(0, 0, 0), mask=""):
    """find closest atom to a given point

    Parameters
    ----------
    top : Topology
    frame : Frame
    point : tuple or 1-D array, default (0, 0, 0)
    mask : str, default ''
        if given, only find closest atom in this mask

    Returns
    -------
    dict : {"index": atom_index, "distance": distance}

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> top = traj.top
    >>> frame = traj[0]
    >>> closest_atom(top=top, frame=frame, point=(0, 0, 0), mask='')
    {'index': 1689, 'distance': 22.492717051651973}
    """
    if mask == '':
        indices = np.arange(0, top.n_atoms)
    else:
        indices = top.select(mask)

    min_dist = 1E6
    min_index = 0
    point = np.asarray(point)

    for index in indices:
        dist = pair_distance(frame.xyz[index], point)
        if dist < min_dist:
            min_dist = dist
            min_index = index

    return dict(index=min_index, distance=min_dist)