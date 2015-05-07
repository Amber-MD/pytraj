"""
pytraj
"""
from __future__ import absolute_import
import os
#from pytraj import base
from .utils.check_and_assert import _import
# create adict and analdict objects here to we can use below in some modules
# >>> from pytraj import adict, analdict 

from .core import Atom, Residue, Molecule
from pytraj.Topology import Topology
from .ArgList import ArgList
#from .NameType import NameType
from .AtomMask import AtomMask
from .math import Vec3
from .CpptrajState import CpptrajState
from .Frame import Frame
from .Trajectory import Trajectory
from .TrajectoryIterator import TrajectoryIterator
from .trajs.Trajout import Trajout
from .datasets.cast_dataset import cast_dataset
from .parms.ParmFile import ParmFile
from . import io

# dataset stuff
from .dataframe import to_dataframe
from .data_sample.load_sample_data import load_sample_data
from .DataSetList import DataSetList
from .DataFileList import DataFileList

# actions and analyses
from .actions import allactions
from .analyses import allanalyses
from ._common_actions import calculate
from . import common_actions
from .action_dict import ActionDict
from .analysis_dict import AnalysisDict
adict = ActionDict()
analdict = AnalysisDict()
from .misc import info
from .run_tests import run_tests

# turn off verbose in cpptraj
# TODO: need to move set_world_silent and set_error_silent to the same file
from ._set_silent import set_error_silent, set_world_silent

set_world_silent(True)

# we still need cpptraj notify error
#set_error_silent(True)
