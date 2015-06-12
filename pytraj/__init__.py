"""
pytraj
"""
from __future__ import absolute_import
import os

from .core import Atom, Residue, Molecule
from pytraj.Topology import Topology
from .ArgList import ArgList
from .AtomMask import AtomMask
from .math import Vec3
from .core.CpptrajState import CpptrajState
from .Frame import Frame
from .Trajectory import Trajectory
from .TrajectoryIterator import TrajectoryIterator
from .trajs.Trajout import Trajout
from .datasets.cast_dataset import cast_dataset
from .parms.ParmFile import ParmFile
from . import io
from .io import  (load, iterload, load_remd, iterload_remd,
                  _load_from_filelist, _iterload_from_filelist,
                  load_pdb_rcsb, load_pdb,
                  load_pseudo_parm, load_cpptraj_file,
                  load_datafile, load_hdf5,
                  load_sample_data,
                  load_ParmEd, load_full_ParmEd,
                  load_mdtraj,
                  load_MDAnalysis, load_MDAnalysisIterator,
                  load_topology, read_parm, write_parm, 
                  save, write_traj,
                  read_pickle, read_json,
                  to_pickle, to_json,
                  )

# dataset stuff
from .data_sample.load_sample_data import load_sample_data
from .DataSetList import DataSetList

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

from ._shared_methods import _frame_iter_master as frame_iter_master

# turn off verbose in cpptraj
# TODO: need to move set_world_silent and set_error_silent to the same file
from ._set_silent import set_error_silent, set_world_silent
from ._set_silent import set_world_silent as set_cpptraj_verbose

set_world_silent(True)
