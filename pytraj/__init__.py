"""
pytraj
"""
from __future__ import absolute_import
from sys import platform as _platform
import sys
import os

from .version import version as __version__
from .c_options import info as compiled_info
from .c_options import __cpptraj_version__
from .c_options import __cpptraj_internal_version__

if 'BINTRAJ' not in compiled_info():
    from warnings import warn
    warn('linking to libcpptraj that were not installed with libnetcdf')

from .c_action.actionlist import ActionList, pipe, compute
Pipeline = ActionList

try:
    import numpy as np
    np.set_printoptions(threshold=10)
except ImportError:
    raise ImportError("require numpy")
    
from . import options

# import partial from functools
from functools import partial

from .externals.six import string_types
from .core import Atom, Residue, Molecule
from .core.c_core import CpptrajState, ArgList, AtomMask
from .core.c_core import Command
from .core.c_core import Command
dispatch = Command.dispatch
from . import array, c_commands
from .topology import Topology, ParmFile
from .math import Vec3
from .frame import Frame
from .shared_methods import iterframe_master
from .trajectory import Trajectory
from .trajectory_iterator import TrajectoryIterator
from .c_traj.c_trajout import TrajectoryWriter
from .datasets.cast_dataset import cast_dataset
from .datasetlist import DatasetList as Dataset
from . import io
from .io import (load,
                 iterload,
                 load_remd,
                 iterload_remd,
                 _load_from_frame_iter,
                 load_pdb_rcsb,
                 load_cpptraj_file,
                 load_sample_data,
                 load_parmed,
                 load_leap,
                 load_antechamber,
                 load_topology,
                 load_batch,
                 write_parm,
                 get_coordinates,
                 save,
                 write_traj,
                 read_pickle,
                 to_pickle,
                 select_atoms,
                 )

from .frameiter import iterframe, iterchunk

# alias
write_trajectory = write_traj

load_from_frame_iter = _load_from_frame_iter

# dataset stuff
from .datafiles import load_cpptraj_state
from .datasetlist import DatasetList

# tool
from . import tools

# actions and analyses
from .c_action import c_action as allactions
from .c_action import c_action
from .c_analysis import c_analysis as allanalyses
from .c_analysis import c_analysis

from .dssp_analysis import calc_dssp, dssp_allatoms, dssp_allresidues
from .nucleic_acid_analysis import nastruct
from .nmr import ired_vector_and_matrix, _ired, NH_order_parameters
from .hbond_analysis import hbond
from .energy_analysis import esander

from .all_actions import (
    calc_rmsd_nofit, rmsd, symmrmsd, rmsd_perres, distance_rmsd, calc_multidihedral,
    autoimage, image, calc_angle, calc_dihedral, calc_distance,
    calc_pairwise_distance, calc_center_of_mass, calc_center_of_geometry,
    calc_jcoupling, calc_surf, calc_molsurf, calc_radgyr, calc_rdf, calc_vector,
    calc_pairwise_rmsd, calc_atomicfluct, atomicfluct, calc_bfactors,
    calc_rotation_matrix, calc_watershell, calc_volume, calc_mindist,
    calc_atomiccorr,
    # lifetime,
    pucker,
    get_average_frame, get_velocity, _dihedral_res,
    native_contacts, principal_axes,
    align_principal_axis,
    timecorr, center, translate, rotate,
    rotate_dihedral, make_structure, scale, randomize_ions,
    set_dihedral, crank, closest, search_neighbors, replicate_cell,
    calc_pairdist, _grid, grid, transform, lowestcurve, calc_diffusion, calc_volmap,
    calc_multivector, pca, projection,
    xcorr, acorr, velocityautocorr,
    check_structure,
    calc_matrix,
    superpose, align, strip,
    density, gist,
    wavelet)

from .matrix import dist as distance_matrix
from . import cluster

from .dihedral_analysis import (calc_phi, calc_psi, calc_alpha, calc_beta,
                        calc_omega, calc_chin, calc_chip, calc_delta,
                        calc_epsilon, calc_gamma, calc_zeta,
                        calc_omega, calc_nu1, calc_nu2)

from .c_action.c_action import ActionDict
from .c_analysis.analysis_dict import AnalysisDict
from . import matrix
from . import dihedral_analysis
from . import vector

# others
from .misc import info
from .run_tests import run_tests


# turn off verbose in cpptraj
# TODO: need to move set_world_silent and set_error_silent to the same file
from .c_options import set_error_silent, set_world_silent
from .cyutils import _fast_iterptr as iterframe_from_array

# create alias
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

# compat with cpptraj, (FIXME)
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

def set_cpptraj_verbose(cm=True):
    if cm:
        set_world_silent(False)
    else:
        set_world_silent(True)

set_world_silent(True)
_verbose = set_cpptraj_verbose

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
           'Trajectory', 'TrajectoryIterator',
           'ActionList', 'ActionDict', 'AnalysisDict', 'adict', 'analdict',
           'dispatch', 'iterchunk', 'iterframe',
           'select', 'set_cpptraj_verbose', 'show_versions',
           'dihedral_analysis', 'hbond_analysis', 'dssp_analysis',
           'nucleic_acid_analysis',])
