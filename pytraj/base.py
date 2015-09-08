"""
import baseclasses for pytraj
"""
from __future__ import absolute_import
# TODO : make this file shorter
from .datasets.cast_dataset import cast_dataset
from .Frame import Frame
from .core.Atom import Atom
from .core.TopologyList import TopologyList
from .datafiles.datafiles import DataFileList
from .core.ActionList import ActionList
from .core.cpptraj_core import CpptrajState
from .datasets.DataSetList import DataSetList

from .core.cpptraj_core import AtomMask
from .api import Trajectory
from .Topology import Topology
from .core.cpptraj_core import ArgList
from .TrajectoryIterator import TrajectoryIterator
from .trajs.Trajout import Trajout
from . import cpptraj_dict

# `Trajectory` is alias of `Trajectory`
__all__ = ['Atom', 'Topology', 'TopologyList', 'Frame', 'Trajectory',
           'AtomMask', 'ArgList', 'CpptrajState', 'TrajectoryIterator',
           'DataSetList', 'DataFileList', 'ActionList', 'Trajout',
           'cast_dataset', 'cpptraj_dict']
