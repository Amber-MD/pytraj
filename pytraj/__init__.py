"""
pytraj
"""
from __future__ import absolute_import
from sys import platform as _platform
import sys
import os

from .version import version as __version__

# checking cpptraj version first
from .cpp_options import info as compiled_info
from .cpp_options import __cpptraj_version__
from .cpp_options import __cpptraj_internal_version__

_v = __cpptraj_internal_version__
# TODO: follow python's rule
if _v < 'V4.2.7':
    raise RuntimeError("need to have cpptraj version >= v4.2.7")

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
    if 'linux' in _platform and not _find_lib("cpptraj"):
        raise ImportError(
            "can not find libcpptraj. Make sure to install it "
            "or export LD_LIBRARY_PATH correctly")

try:
    from .core import Atom, Residue, Molecule
    from .c_action.actionlist import ActionList, pipe, do
    compute = do
    Pipeline = ActionList

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

from . import options

# import partial from functools
from functools import partial

from .externals.six import string_types
from .core import Atom, Residue, Molecule
from .core.c_core import CpptrajState, ArgList, AtomMask, _load_batch
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
                 load_ParmEd,
                 load_topology,
                 read_parm,
                 write_parm,
                 get_coordinates,
                 save,
                 write_traj,
                 read_pickle,
                 read_json,
                 to_pickle,
                 to_json, )

# alias
write_trajectory = write_traj

load_from_frame_iter = _load_from_frame_iter

# dataset stuff
from .datafiles import load_cpptraj_state
from .datasetlist import DatasetList

# alias
load_cpptrajstate = load_cpptraj_state
load_state = load_cpptraj_state

# tool
from . import tools

# actions and analyses
from .c_action import c_action as allactions
from .c_action import c_action
from .c_analysis import c_analysis as allanalyses
from .c_analysis import c_analysis
from . import all_actions

from .dssp_ import calc_dssp, dssp_allatoms, dssp_allresidues
from .nucleic_acid_ import nastruct
from .nmr import ired_vector_and_matrix, _ired, NH_order_parameters
from .hbond_ import hbond
from .externals.energy import energy_decomposition

from .all_actions import (
    calc_rmsd_nofit, rmsd, rmsd_perres, distance_rmsd, calc_multidihedral,
    autoimage, calc_angle, calc_dihedral, calc_distance,
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
    rotate_dihedral, make_structure, scale, clustering_dataset, randomize_ions,
    set_dihedral, crank, closest, search_neighbors, replicate_cell,
    calc_pairdist, _grid, transform, lowestcurve, calc_diffusion, calc_volmap,
    calc_multivector, pca, projection,
    xcorr, acorr,
    check_structure,
    calc_matrix)

from .matrix import dist
distance_matrix = dist
from . import cluster

from .dihedral_ import (calc_phi, calc_psi, calc_alpha, calc_beta,
                        calc_omega, calc_chin, calc_chip, calc_delta,
                        calc_epsilon, calc_gamma, calc_zeta,
                        calc_omega, calc_nu1, calc_nu2)

from .c_action.action_dict import ActionDict
from .c_analysis.analysis_dict import AnalysisDict
from . import matrix
from . import dihedral_
from . import vector

# others
from .misc import info
from .run_tests import run_tests


# turn off verbose in cpptraj
# TODO: need to move set_world_silent and set_error_silent to the same file
from .cpp_options import set_error_silent, set_world_silent
from .cyutils import _fast_iterptr as iterframe_from_array

# create alias
check_overlap = check_structure
checkoverlap = check_structure
fetch_pdb = load_pdb_rcsb
rmsd_nofit = calc_rmsd_nofit
drmsd = distance_rmsd
search_hbonds = hbond
distance = calc_distance
distances = calc_distance
pairwise_distance = calc_pairwise_distance
angle = calc_angle
angles = calc_angle
dihedral = calc_dihedral
dihedrals = calc_dihedral
jcoupling = calc_jcoupling
rmsf = calc_atomicfluct
pairwise_rmsd = calc_pairwise_rmsd
rms2d = calc_pairwise_rmsd
rotation_matrix = calc_rotation_matrix
multidihedral = calc_multidihedral
dssp = calc_dssp
bfactors = calc_bfactors
volume = calc_volume
radgyr = calc_radgyr
rdf = calc_rdf
pairdist = calc_pairdist
multivector = calc_multivector
atomiccorr = calc_atomiccorr
surf = calc_surf
molsurf = calc_molsurf
center_of_mass = calc_center_of_mass
center_of_geometry = calc_center_of_geometry
watershell = calc_watershell
mean_structure = get_average_frame
average_frame = get_average_frame
load_parmed = load_ParmEd
mindist = calc_mindist
# compat with cpptraj
nativecontacts = native_contacts
pair_distribution = pairdist
lowest_curve = lowestcurve
diffusion = calc_diffusion
volmap = calc_volmap
calc_pca = pca
randomizeions = randomize_ions

adict = ActionDict()
analdict = AnalysisDict()

# import parallel package after all pytraj or cpptraj's method so we
# can import them to parallel namespace
# import _pmap here to be called from nmr module
from .parallel.multiprocessing_ import pmap, _pmap
from .parallel.mpi import pmap_mpi
from .parallel.base import _load_batch_pmap
from . import view


def load_batch(traj, txt):
    '''perform calculation for traj with cpptraj's batch style. This is for internal use.

    Parameters
    ----------
    traj : pytraj.TrajectoryIterator
    txt : text or a list of test
        cpptraj's commands

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> text = """
    ... autoimage
    ... radgyr @CA nomax
    ... molsurf !@H=
    ... """
    >>> state = pt.load_batch(traj, text)
    >>> state = state.run()
    >>> state.data
    <pytraj.datasets.CpptrajDatasetList - 3 datasets>

    >>> # raise if not TrajectoryIterator
    >>> traj2 = pt.Trajectory(xyz=traj.xyz, top=traj.top)
    >>> not isinstance(traj2, pt.TrajectoryIterator)
    True
    >>> pt.load_batch(traj2, text)
    Traceback (most recent call last):
        ...
    ValueError: only support TrajectoryIterator
    '''
    if not isinstance(traj, TrajectoryIterator):
        raise ValueError('only support TrajectoryIterator')
    return _load_batch(txt, traj=traj)


load_pipeline = load_batch


def superpose(traj, *args, **kwd):
    '''

    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_ala3()[:]
    >>> traj = pt.superpose(traj)
    >>> isinstance(traj, pt.Trajectory)
    True
    '''
    traj.superpose(*args, **kwd)
    return traj


def set_cpptraj_verbose(cm=True):
    if cm:
        set_world_silent(False)
    else:
        set_world_silent(True)


set_world_silent(True)
_verbose = set_cpptraj_verbose


def iterframe(traj, *args, **kwd):
    """create frame iterator with given indices, mask or some iter_options

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> for frame in pt.iterframe(traj, 0, 8, 2): pass
    >>> for frame in pt.iterframe(traj, 4, mask='@CA'): print(frame)
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>

    # create frame iterator for given indices
    >>> for frame in pt.iterframe(traj, frame_indices=[0, 7, 3]): print(frame)
    <Frame with 5293 atoms>
    <Frame with 5293 atoms>
    <Frame with 5293 atoms>


    >>> fi = pt.iterframe(traj)
    >>> # iterframe its self
    >>> fi = pt.iterframe(fi)

    See also
    --------
    pytraj.TrajectoryIterator.iterframe
    """
    if hasattr(traj, 'iterframe'):
        return traj.iterframe(*args, **kwd)
    else:
        return iterframe_master(traj)


def iterchunk(traj, *args, **kwd):
    """iterate ``traj`` by chunk

    Parameters
    ----------
    traj : TrajectoryIterator
    chunksize : int
        the number of frames in each chunk
    start : int, default 0
        start frame to iterate
    start : int, default -1 (last frame)
        stop frame
    autoimage : bool, default False
        if True, do autoimage for chunk

    Return
    ------
    pytraj.Trajectory, n_frames=chunksize
        The final chunk might not have the n_frames=chunksize

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> for frame in pt.iterchunk(traj, 4): pass
    >>> for frame in pt.iterchunk(traj, chunksize=4, start=2): pass
    >>> for frame in pt.iterchunk(traj, chunksize=4, start=2, stop=9): pass
    >>> for frame in pt.iterchunk(traj, chunksize=4, autoimage=True): pass

    See also
    --------
    pytraj.TrajectoryIterator.iterchunk
    """
    return traj.iterchunk(*args, **kwd)


def select_atoms(mask, topology):
    '''return atom indices

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> atom_indices = pt.select_atoms('@CA', traj.top)
    >>> atom_indices
    array([  4,  15,  39, ..., 159, 173, 197])
    >>> pt.select_atoms(traj.top, '@CA')
    array([  4,  15,  39, ..., 159, 173, 197])
    '''
    if isinstance(mask, Topology) and isinstance(topology, string_types):
        mask, topology = topology, mask
    return topology.select(mask)


select = select_atoms


def strip(obj, mask):
    '''return a new Trajectory or Topology with given mask.

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> traj.n_atoms
    5293
    >>> pt.strip(traj, '!@CA').n_atoms
    12
    >>> pt.strip(traj.top, '!@CA').n_atoms
    12
    >>> pt.strip('!@CA', traj.top).n_atoms
    12
    >>> t0 = traj[:3]
    >>> pt.strip(t0, '!@CA').n_atoms
    12
    >>> fi = pt.iterframe(traj, stop=3)
    >>> fi = pt.strip(fi, '!@CA')
    >>> for frame in fi:
    ...    print(frame)
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>

    >>> # raise ValueError
    >>> pt.strip(0, '@CA')
    Traceback (most recent call last):
        ...
    ValueError: object must be either Trajectory or Topology
    '''

    if isinstance(obj, string_types) and not isinstance(
            mask, string_types):
        obj, mask = mask, obj

    kept_mask = '!(' + mask + ')'

    if isinstance(obj, (Topology, Trajectory)):
        # return new Topology or new Trajectory
        return obj[kept_mask]
    elif isinstance(obj, TrajectoryIterator):
        # return a FrameIterator
        return obj(mask=kept_mask)
    elif hasattr(obj, 'mask'):
        obj.mask = kept_mask
        return obj
    else:
        raise ValueError('object must be either Trajectory or Topology')


def run(fi):
    '''shortcut for `for frame in fi: pass`

    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> fi = pt.pipe(traj, ['autoimage', 'rms', 'center :1-13'])
    >>> pt.run(fi)
    '''
    for _ in fi:
        pass


def show():
    # just delay importing
    """show plot
    """
    from matplotlib import pyplot
    pyplot.show()


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
