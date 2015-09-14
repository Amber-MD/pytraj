"""
pytraj
"""
from __future__ import absolute_import
from sys import platform as _platform
import os

# checking cpptraj version first
from .cpp_options import info as compiled_info
from .cpp_options import __cpptraj_version__
from .cpp_options import __cpptraj_internal_version__

_v = __cpptraj_internal_version__
# TODO: follow python's rule
if _v < 'V4.2':
    raise RuntimeError("need to have cpptraj version >= v4.2")

if 'BINTRAJ' not in compiled_info():
    from warnings import warn
    warn('linking to libcpptraj that were not installed with libnetcdf')

from .tools import find_lib as _find_lib

# check `libcpptraj` and raise ImportError
# only check for Linux since I don't know much about
# OS X and Windows
try:
    # try to check `libcpptraj` that not in LD_LIBRARY_PATH search
    # in _find_lib
    from .core import Atom
except ImportError:
    if 'linux' in _platform:
        if not _find_lib("cpptraj"):
            raise ImportError(
                "can not find libcpptraj. Make sure to install it "
                "or export LD_LIBRARY_PATH correctly")

try:
    from .core import Atom, Residue, Molecule
except ImportError:
    import os
    source_folders = ['./scripts', './devtools', './docs']
    is_source_folder = True
    for f in source_folders:
        is_source_folder = False if not os.path.exists(f) else True
    if is_source_folder:
        raise ImportError("you are import pytraj in source folder. "
                          "Should go to another location and try again")
try:
    import numpy as np
    np.set_printoptions(threshold=10)
except ImportError:
    np = None

from .__version__ import __version__
version = __version__
from . import options

# import partial from functools
from functools import partial

from .core import Atom, Residue, Molecule
from .core.cpptraj_core import CpptrajState, ArgList, AtomMask
from . import array
from .Topology import Topology, ParmFile
from .math import Vec3
from .Frame import Frame
from .api import Trajectory
from .TrajectoryIterator import TrajectoryIterator
from .trajs.Trajout import Trajout
from .datasets.cast_dataset import cast_dataset
from .datasetlist import DatasetList as Dataset
from . import io
from .io import (load, iterload, load_remd, iterload_remd, _load_from_frame_iter, load_pdb_rcsb,
                 load_pdb, load_cpptraj_file, load_sample_data,
                 load_ParmEd, load_mdtraj, load_MDAnalysis,
                 load_MDAnalysisIterator, load_topology, read_parm, write_parm,
                 get_coordinates, save, write_traj, read_pickle, read_json,
                 to_pickle, to_json, )

load_from_frame_iter = _load_from_frame_iter

# dataset stuff
from .datafiles.load_sample_data import load_sample_data
from .datasetlist import DatasetList

# tool
from . import tools

# actions and analyses
from .actions import CpptrajActions as allactions
from .analyses import CpptrajAnalyses as allanalyses
from . import common_actions
from .dssp_analysis import calc_dssp
from .common_actions import (
    rmsd, rmsd_perres, search_hbonds,
    calc_multidihedral, autoimage, nastruct, calc_angle, calc_dihedral,
    calc_distance, calc_center_of_mass, calc_center_of_geometry, calc_dssp,
    calc_jcoupling, calc_molsurf, calc_radgyr, calc_rdf, calc_vector,
    calc_pairwise_rmsd, calc_atomicfluct, calc_bfactors, calc_density,
    calc_watershell, calc_volume, calc_mindist, lifetime, get_average_frame,
    _dihedral_res, energy_decomposition, native_contacts,
    auto_correlation_function, principal_axes, cross_correlation_function,
    timecorr, center, translate, rotate, rotate_dihedral, make_structure,
    scale, do_clustering, clustering_dataset, _rotate_dih, randomize_ions,
    crank, closest, search_neighbors, replicate_cell)

# create alias
distance = calc_distance
angle = calc_angle
dihedral = calc_dihedral
jcoupling = calc_jcoupling
nucleic_acid_analysis = nastruct
calc_RMSF = calc_atomicfluct
pairwise_rmsd = calc_pairwise_rmsd
multidihedral = calc_multidihedral
xcorr = cross_correlation_function
acorr = auto_correlation_function
dssp = calc_dssp
bfactors = calc_bfactors
density = calc_density
volume = calc_volume
radgyr = calc_radgyr
molsurf = calc_molsurf
center_of_mass = calc_center_of_mass
center_of_geometry = calc_center_of_geometry
watershell = calc_watershell
get_average_structure = get_average_frame
average_frame = get_average_frame
load_parmed = load_ParmEd
from_parmed = load_ParmEd
clustering = do_clustering
mindist = calc_mindist

from .matrix_analysis import distance_matrix
from .dihedral_analysis import (
    calc_phi, calc_psi, calc_alpha, calc_beta, calc_omega, calc_chin,
    calc_chip, calc_delta, calc_epsilon, calc_gamma, calc_zeta, calc_omega,
    calc_nu1, calc_nu2)

from .action_dict import ActionDict
from .analysis_dict import AnalysisDict
adict = ActionDict()
analdict = AnalysisDict()
from . import matrix_analysis
from . import dihedral_analysis
from . import vector_analysis

# others
from .misc import info
from .run_tests import run_tests

# parallel
from .pmap import pmap

from ._shared_methods import iterframe_master

# turn off verbose in cpptraj
# TODO: need to move set_world_silent and set_error_silent to the same file
from .cpp_options import set_error_silent, set_world_silent


def to_numpy_Trajectory(traj, top, unitcells=None):
    # TODO: move to `io`?
    from . import api
    import numpy as np
    from ._xyz import XYZ

    t = api.Trajectory(top=top)
    if isinstance(traj, np.ndarray) or isinstance(traj, XYZ):
        t.xyz = np.asarray(traj)
    elif hasattr(traj, 'xyz'):
        t.xyz = traj.xyz
        if hasattr(traj, 'unitcells'):
            t.unitcells = traj.unitcells
    else:
        t.xyz = get_coordinates(traj)
    if unitcells is not None:
        t.unitcells = unitcells
    return t


def to_mdtraj(traj, top=None):
    # TODO: move to `io`?
    from pytraj.utils.context import goto_temp_folder
    import mdtraj as md

    _top = top if top is not None else traj.top
    xyz = get_coordinates(traj)

    with goto_temp_folder():
        _top.save("tmp.prmtop")
        top = md.load_prmtop("tmp.prmtop")
        return md.Trajectory(xyz, top)


def set_cpptraj_verbose(cm=True):
    if cm:
        set_world_silent(False)
    else:
        set_world_silent(True)


set_world_silent(True)
_verbose = set_cpptraj_verbose


def iterframe(traj, *args, **kwd):
    """

    Examples
    --------
    >>> import pytraj as pt
    >>> for frame in pt.iterframe(traj, 0, 8, 2): print(frame)
    >>> for frame in pt.iterframe(traj, 0, 8, 2, mask='@CA'): print(frame)
    """
    return traj.iterframe(*args, **kwd)

from ._cyutils import _fast_iterptr as iterframe_from_array


def iterchunk(traj, *args, **kwd):
    """

    Examples
    --------
    >>> import pytraj as pt
    >>> for chunk in pt.iterchunk(traj, 4): print(chunk)
    >>> for chunk in pt.iterframe(traj, 4, mask='@CA'): print(chunk)
    """
    return traj.iterchunk(*args, **kwd)


def show():
    # just delay importing
    """show plot
    """
    from matplotlib import pyplot
    pyplot.show()


def savefig(fname, *args, **kwd):
    from matplotlib import pyplot
    pyplot.savefig(fname, *args, **kwd)


def show_versions():
    """
    """
    print("pytraj version = ", version)
    print("cpptraj version = ", __cpptraj_version__)
    print("cpptraj internal version = ", __cpptraj_internal_version__)
    print("cpptraj compiled flag = ", compiled_info())
