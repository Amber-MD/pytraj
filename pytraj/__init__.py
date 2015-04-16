"""
pytraj
"""
from __future__ import absolute_import
import os
#from pytraj import base
from .utils.check_and_assert import _import

from .import io

from .action_dict import ActionDict
from .analysis_dict import AnalysisDict
# create adict and analdict objects here to we can use below in some modules
# >>> from pytraj import adict, analdict 
adict = ActionDict()
analdict = AnalysisDict()

from .run_tests import run_tests
from .Atom import Atom
from .Frame import Frame
from .FrameArray import FrameArray
from .Topology import Topology
from .NameType import NameType
from .ArgList import ArgList
from .AtomMask import AtomMask
from .CpptrajState import CpptrajState
from .TrajReadOnly import TrajReadOnly
from .trajs.Trajout import Trajout
from .cast_dataset import cast_dataset
from .parms.ParmFile import ParmFile
from .misc import info

# dataset stuff
from .dataframe import to_dataframe
from .data_sample.load_sample_data import load_sample_data
from .DataSetList import DataSetList
from .DataFileList import DataFileList

# actions
from .actions import allactions
from .analyses import allanalyses
from ._common_actions import calculate
from . import common_actions

# turn off verbose in cpptraj
# TODO: need to move set_world_silent and set_error_silent to the same file
from ._utils import set_world_silent
from ._set_silent import set_error_silent

set_world_silent(True)

# we still need cpptraj notify error
#set_error_silent(True)
