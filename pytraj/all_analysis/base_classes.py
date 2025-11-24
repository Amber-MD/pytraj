"""
Base classes and utilities used across all analysis modules
"""
from __future__ import absolute_import
import numpy as np
try:
    from enum import StrEnum
except ImportError:
    from enum import Enum
    class StrEnum(str, Enum):
        def __new__(cls, value):
            member = str.__new__(cls, value)
            member._value_ = value
            return member

from typing import Any, Callable, List, Union
from functools import partial

from ..utils.get_common_objects import (
    get_topology,
    get_data_from_dtype,
    get_list_of_commands,
    get_reference,
    get_fiterator,
    super_dispatch,
)
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
from ..core.c_core import CpptrajState, Command


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


class AnalysisRunner:
    """Class to handle analysis calculations with cpptraj"""

    def __init__(self):
        self.state = CpptrajState()
        self.dslist = CpptrajDatasetList()

    def run(self, command, traj, top=None):
        """Run an analysis command on trajectory data"""
        if top is not None:
            self.state.load(top)

        # Process trajectory frames
        for frame in traj:
            # Run analysis on frame
            pass

        return self.dslist


class CommandBuilder:
    """Helper class to build cpptraj command strings"""

    def __init__(self):
        self.commands = []

    def add_command(self, action, mask="", options=""):
        """Add a command to the command list"""
        cmd_parts = [action]
        if mask:
            cmd_parts.append(mask)
        if options:
            cmd_parts.append(options)

        command = " ".join(cmd_parts)
        self.commands.append(command)
        return self

    def build(self):
        """Build the final command string"""
        return "\n".join(self.commands)

    def clear(self):
        """Clear all commands"""
        self.commands.clear()


class CommandType(StrEnum):
    INT = 'int'
    STR = 'str'
    LIST = 'list'


def _check_command_type(command):
    """Check the type of command and return appropriate CommandType"""
    if isinstance(command, int):
        return CommandType.INT
    elif isinstance(command, str):
        return CommandType.STR
    elif isinstance(command, (list, tuple)):
        return CommandType.LIST
    else:
        raise ValueError(f"Unsupported command type: {type(command)}")


def _assert_mutable(trajiter):
    """make sure the input is not TrajectoryIterator"""
    if isinstance(trajiter, TrajectoryIterator):
        raise ValueError("TrajectoryIterator is immutable")


def _calculate_angles_for_int_array(traj, integer_array, n_frames, dtype):
    """Internal function to calculate angles using integer arrays"""
    if integer_array.ndim == 1:
        integer_array = np.atleast_2d(integer_array)

    arr = np.empty([n_frames, len(integer_array)])

    for idx, frame in enumerate(iterframe_master(traj)):
        arr[idx] = frame._angle(integer_array)

    arr = arr.T
    if dtype == 'ndarray':
        return arr
    else:
        dslist = DatasetList({'angle': arr})
        return get_data_from_dtype(dslist, dtype)


def _create_and_compute_action_list(list_of_commands: List[str],
                                    top,
                                    traj: 'Trajectory',
                                    action: Callable,
                                    dtype: str,
                                    args: tuple,
                                    kwargs: dict) -> Any:
    """Create and compute an action list for multiple commands"""
    cpptraj_action_datasets = CpptrajDatasetList()
    action_list = ActionList()

    for command in list_of_commands:
        action_list.add(action(), command, top, dslist=cpptraj_action_datasets, *args, **kwargs)

    action_list.compute(traj)
    return get_data_from_dtype(cpptraj_action_datasets, dtype)


# Export commonly used functions and classes
__all__ = [
    'DatasetType',
    'AnalysisRunner',
    'CommandBuilder',
    'CommandType',
    '_check_command_type',
    '_assert_mutable',
    '_calculate_angles_for_int_array',
    '_create_and_compute_action_list',
]