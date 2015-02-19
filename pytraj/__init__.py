"""
pytraj
"""
import os
from pytraj import base
from pytraj.utils.check_and_assert import _import

# TODO : should we load those module here or use .base?
from pytraj import io

from pytraj.action_dict import ActionDict
from pytraj.analysis_dict import AnalysisDict
# make adict instance here so run_tests can use it.
adict = ActionDict()
analdict = AnalysisDict()

from pytraj.run_tests import run_tests
from pytraj.AtomSelect import AtomSelect
from pytraj.Frame import Frame
from pytraj.FrameArray import FrameArray
from pytraj.Topology import Topology
from pytraj.NameType import NameType
from pytraj.ArgList import ArgList
from pytraj.AtomMask import AtomMask
from pytraj.CpptrajState import CpptrajState
from pytraj.TrajReadOnly import TrajReadOnly
from pytraj.trajs.Trajout import Trajout
from pytraj.cast_dataset import cast_dataset
from pytraj.parms.ParmFile import ParmFile
from pytraj.misc import action_help
from pytraj.data_sample.load_sample_data import load_sample_data
from pytraj.DataSetList import DataSetList
from pytraj.DataFileList import DataFileList

# actions
from pytraj.actions import allactions
from pytraj.analyses import allanalyses
from pytraj.misc import calculate

#try:
#    amber_home = os.environ['AMBERHOME']
#except:
#    raise EnvironmentError("must set AMBERHOME")

# make adict instance of ActionDict
