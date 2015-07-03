"""
pytraj
"""
from __future__ import absolute_import

from .__version__ import __version__
version = __version__
from . import options

# import partial from functools
from functools import partial

from .core import Atom, Residue, Molecule
from .core.CpptrajState import CpptrajState
from .import array
from .Topology import Topology
from .ArgList import ArgList
from .AtomMask import AtomMask
from .math import Vec3
from .Frame import Frame
from .Trajectory import Trajectory
from .TrajectoryIterator import TrajectoryIterator
from .trajs.Trajout import Trajout
from .datasets.cast_dataset import cast_dataset
from .parms.ParmFile import ParmFile
from . import io
from .io import (load, iterload, load_remd, iterload_remd,
                 _load_from_filelist, _iterload_from_filelist,
                 _load_from_frame_iter,
                 load_pdb_rcsb, load_pdb,
                 load_pseudo_parm, load_cpptraj_file,
                 load_datafile, load_hdf5,
                 load_sample_data,
                 load_ParmEd, load_full_ParmEd,
                 load_mdtraj,
                 load_MDAnalysis, load_MDAnalysisIterator,
                 load_topology, read_parm, write_parm,
                 get_coordinates,
                 save, write_traj,
                 read_pickle, read_json,
                 to_pickle, to_json,
                 )

load_from_frame_iter = _load_from_frame_iter

# dataset stuff
from .datafiles.load_sample_data import load_sample_data
from .datasetlist import DatasetList

# tool
from . import tools

# actions and analyses
from .actions import CpptrajActions as allactions
from .analyses import CpptrajAnalyses as allanalyses
from ._common_actions import calculate
from . import common_actions
from . dssp_analysis import calc_dssp
from . common_actions import (rmsd, search_hbonds,
                              calc_multidihedral,
                              autoimage, nastruct,
                              calc_angle, calc_dihedral, calc_distance,
                              calc_center_of_mass, calc_center_of_geometry,
                              calc_dssp, calc_jcoupling, calc_molsurf,
                              calc_radgyr, calc_rdf, calc_vector,
                              calc_pairwise_rmsd,
                              calc_atomicfluct,
                              calc_bfactors,
                              energy_decomposition,)

# create alias
nucleic_acid_analysis = nastruct
calc_RMSF = calc_atomicfluct

from . matrix_analysis import distance_matrix
from . dihedral_analysis import (
    calc_phi, calc_psi, calc_omega, calc_chin, calc_chip)

from .action_dict import ActionDict
from .analysis_dict import AnalysisDict
adict = ActionDict()
analdict = AnalysisDict()
from . import matrix_analysis
from . import dihedral_analysis

# others
from .misc import info
from .run_tests import run_tests

from ._shared_methods import _frame_iter_master as frame_iter_master

# turn off verbose in cpptraj
# TODO: need to move set_world_silent and set_error_silent to the same file
from ._set_silent import set_error_silent, set_world_silent
from ._set_silent import set_world_silent as set_cpptraj_verbose

set_world_silent(True)


def show_versions():
    """
    """
    from .__cpptraj_version__ import info as compiled_info
    from .__cpptraj_version__ import __cpptraj_version__
    from .__cpptraj_version__ import __cpptraj_internal_version__
    print("pytraj version = ", version)
    print("cpptraj version = ", __cpptraj_version__)
    print("cpptraj internal version = ", __cpptraj_internal_version__)
    print("cpptraj compiled flag = ", compiled_info())
