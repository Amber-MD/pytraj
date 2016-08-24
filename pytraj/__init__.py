"""
pytraj
"""
from __future__ import absolute_import

try:
    import numpy as np
    np.set_printoptions(threshold=10)
except ImportError:
    raise ImportError("require numpy")

import sys
import os

from .version import version as __version__
from .core.c_options import info as compiled_info
from .core.c_options import __cpptraj_version__
from .core.c_options import __cpptraj_internal_version__
from .analysis.c_action.actionlist import ActionList
from .analysis.c_action.actionlist import ActionList as Pipeline
from .analysis.c_action.actionlist import pipe
from .analysis.c_action.actionlist import compute
from .core import Atom
from .core import Residue
from .core import Molecule
from .core.c_core import CpptrajState
from .core.c_core import ArgList
from .core.c_core import AtomMask
from .core.c_core import Command
from .datasets import array
from .utils import c_commands
from .utils import tools
from .trajectory.trajectory import Trajectory
from .trajectory.trajectory_iterator import TrajectoryIterator
from .trajectory.c_traj.c_trajout import TrajectoryWriter
from .trajectory.frame import Frame
from .trajectory import frame
from .topology.topology import Topology
from .topology.topology import ParmFile
from .math import Vec3
from .trajectory.shared_methods import iterframe_master
from .trajectory.frameiter import iterframe
from .trajectory.frameiter import iterchunk
from .trajectory.frameiter import FrameIterator
from .datasets.cast_dataset import cast_dataset
from .datasets.datasetlist import DatasetList as Dataset

from . import io
from .io import load
from .io import iterload
from .io import load_remd
from .io import iterload_remd
from .io import _load_from_frame_iter as load_from_frame_iter
from .io import load_pdb_rcsb
from .io import load_cpptraj_file
from .io import load_sample_data
from .io import load_parmed
from .io import load_leap
from .io import load_antechamber
from .io import load_topology
from .io import load_batch
from .io import write_parm
from .io import get_coordinates
from .io import save
from .io import write_traj
from .io import read_pickle
from .io import to_pickle
from .io import select_atoms

# dataset stuff
from .datafiles import load_cpptraj_state
from .datasets.datasetlist import DatasetList

# tool
from .utils import tools

# actions and analyses
from .analysis.c_action import c_action as allactions
from .analysis.c_action import c_action
from .analysis.c_analysis import c_analysis as allanalyses
from .analysis.c_analysis import c_analysis

from .analysis.dssp_analysis import calc_dssp
from .analysis.dssp_analysis import dssp_allatoms
from .analysis.dssp_analysis import dssp_allresidues
from .analysis.energy_analysis import esander
from .analysis.hbond_analysis import hbond
from .analysis.nucleic_acid_analysis import nastruct
from .analysis.nmr import ired_vector_and_matrix
from .analysis.nmr import _ired
from .analysis.nmr import NH_order_parameters

from .analysis import dssp_analysis
from .analysis import energy_analysis
from .analysis import hbond_analysis
from .analysis import nucleic_acid_analysis

from .all_actions import acorr
from .all_actions import align
from .all_actions import align_principal_axis
from .all_actions import atomicfluct
from .all_actions import autoimage
from .all_actions import calc_angle
from .all_actions import calc_atomiccorr
from .all_actions import calc_atomicfluct
from .all_actions import calc_bfactors
from .all_actions import calc_center_of_geometry
from .all_actions import calc_center_of_mass
from .all_actions import calc_diffusion
from .all_actions import calc_dihedral
from .all_actions import calc_distance
from .all_actions import calc_jcoupling
from .all_actions import calc_matrix
from .all_actions import calc_mindist
from .all_actions import calc_molsurf
from .all_actions import calc_multidihedral
from .all_actions import calc_multivector
from .all_actions import calc_pairdist
from .all_actions import calc_pairwise_distance
from .all_actions import calc_pairwise_rmsd
from .all_actions import calc_radgyr
from .all_actions import calc_rdf
from .all_actions import calc_rmsd_nofit
from .all_actions import calc_rotation_matrix
from .all_actions import calc_surf
from .all_actions import calc_vector
from .all_actions import calc_volmap
from .all_actions import calc_volume
from .all_actions import calc_watershell
from .all_actions import center
from .all_actions import check_structure
from .all_actions import closest
from .all_actions import crank
from .all_actions import density
from .all_actions import _dihedral_res
from .all_actions import distance_rmsd
from .all_actions import get_average_frame
from .all_actions import get_velocity
from .all_actions import gist
from .all_actions import grid
from .all_actions import _grid
from .all_actions import image
# from .all_actions import lifetime
from .all_actions import lowestcurve
from .all_actions import make_structure
from .all_actions import native_contacts
from .all_actions import pca
from .all_actions import principal_axes
from .all_actions import projection
from .all_actions import pucker
from .all_actions import randomize_ions
from .all_actions import replicate_cell
from .all_actions import rmsd
from .all_actions import rmsd_perres
from .all_actions import rotate
from .all_actions import rotate_dihedral
from .all_actions import scale
from .all_actions import search_neighbors
from .all_actions import set_dihedral
from .all_actions import strip
from .all_actions import superpose
from .all_actions import symmrmsd
from .all_actions import timecorr
from .all_actions import transform
from .all_actions import translate
from .all_actions import velocityautocorr
from .all_actions import wavelet
from .all_actions import xcorr

from .analysis.matrix import dist as distance_matrix
from .analysis import matrix
from .analysis import vector
from . import cluster

from .analysis import dihedral_analysis
from .analysis.dihedral_analysis import calc_alpha
from .analysis.dihedral_analysis import calc_beta
from .analysis.dihedral_analysis import calc_chin
from .analysis.dihedral_analysis import calc_chip
from .analysis.dihedral_analysis import calc_delta
from .analysis.dihedral_analysis import calc_epsilon
from .analysis.dihedral_analysis import calc_gamma
from .analysis.dihedral_analysis import calc_nu1
from .analysis.dihedral_analysis import calc_nu2
from .analysis.dihedral_analysis import calc_omega
from .analysis.dihedral_analysis import calc_omega
from .analysis.dihedral_analysis import calc_phi
from .analysis.dihedral_analysis import calc_psi
from .analysis.dihedral_analysis import calc_zeta

from .analysis.c_action.c_action import ActionDict
from .analysis.c_analysis.analysis_dict import AnalysisDict

# others
from .utils.misc import info
from .testing.run_tests import run_tests

# turn off verbose in cpptraj
from .core.c_options import set_error_silent
from .core.c_options import set_world_silent
from .core.c_options import set_cpptraj_verbose
from .core.c_options import set_cpptraj_verbose as _verbose
set_world_silent(True)

from .utils.cyutils import _fast_iterptr as iterframe_from_array

# alias
write_trajectory = write_traj
select = select_atoms
dispatch = Command.dispatch
energy_decomposition = esander
check_overlap = check_structure
fetch_pdb = load_pdb_rcsb
rmsd_nofit = calc_rmsd_nofit
search_hbonds = hbond
distances = calc_distance
pairwise_distance = calc_pairwise_distance
angles = calc_angle
dihedrals = calc_dihedral
rmsf = calc_atomicfluct
pairwise_rmsd = calc_pairwise_rmsd
rotation_matrix = calc_rotation_matrix
multidihedral = calc_multidihedral
bfactors = calc_bfactors
rdf = calc_rdf
atomiccorr = calc_atomiccorr
center_of_mass = calc_center_of_mass
center_of_geometry = calc_center_of_geometry
mean_structure = get_average_frame
average_frame = get_average_frame
calc_pca = pca
pair_distribution = pairdist = calc_pairdist

# compat with cpptraj
distance = calc_distance
angle = calc_angle
dihedral = calc_dihedral
jcoupling = calc_jcoupling
dssp = calc_dssp
drmsd = distance_rmsd
checkoverlap = check_structure
radgyr = calc_radgyr
nativecontacts = native_contacts
mindist = calc_mindist
lowest_curve = lowestcurve
diffusion = calc_diffusion
multivector = calc_multivector
volmap = calc_volmap
randomizeions = randomize_ions
molsurf = calc_molsurf
surf = calc_surf
watershell = calc_watershell
pairdist = calc_pairdist
volume = calc_volume
rms2d = calc_pairwise_rmsd

adict = ActionDict()
analdict = AnalysisDict()

# import parallel package after all pytraj or cpptraj's method so we
# can import them to parallel namespace
# import _pmap here to be called from nmr module
from .parallel.multiprocess import pmap, _pmap
from .parallel.mpi import pmap_mpi
from .parallel.base import _load_batch_pmap
from .visualization import view

def show_versions():
    """
    >>> show_versions() # doctest: +SKIP
    """
    print(sys.version)
    print('')
    print("pytraj version = ", __version__)
    print("cpptraj version = ", __cpptraj_version__)
    print("cpptraj internal version = ", __cpptraj_internal_version__)
    print("cpptraj compiled flag = ", compiled_info())


# for website
# do not put __all__ in the top of this file to avoid circular import (all_actions)
__all__ = (io.__all__ 
        + all_actions.__all__
        + dihedral_analysis.__all__
        + ['nastruct']
        + ['esander']
        + ['Atom', 'Residue', 'Molecule', 'Topology', 'Frame', 'AtomMask',
           'Trajectory', 'TrajectoryIterator', 'TrajectoryWriter',
           'ActionList', 'ActionDict', 'AnalysisDict', 'adict', 'analdict',
           'dispatch', 'iterchunk', 'iterframe',
           'select', 'set_cpptraj_verbose', 'show_versions',
           'dihedral_analysis', 'hbond_analysis', 'dssp_analysis',
           'nucleic_acid_analysis','tools'])
