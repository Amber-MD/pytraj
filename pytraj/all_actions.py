from __future__ import absolute_import
import numpy as np
try:
    from enum import StrEnum
except ImportError:
    from enum import Enum
    class StrEnum(str, Enum):
        """
        Enum where members are also (and must be) strs.
        """
        def __new__(cls, value):
            member = str.__new__(cls, value)
            member._value_ = value
            return member

from typing import Any, Callable, List, Union
from functools import partial

from .utils.get_common_objects import (
    get_topology,
    get_data_from_dtype,
    get_list_of_commands,
    get_reference,
    get_fiterator,
    super_dispatch, )
from .utils import ensure_not_none_or_string
from .utils import is_int
from .utils.context import tempfolder
from .utils.context import capture_stdout
from .utils.convert import array_to_cpptraj_atommask
from .utils.convert import array2d_to_cpptraj_maskgroup
from .datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from .datasets.datasetlist import DatasetList
from .trajectory.shared_methods import iterframe_master
from .trajectory.frame import Frame
from .trajectory.trajectory import Trajectory
from .trajectory.trajectory_iterator import TrajectoryIterator
from .utils.decorators import register_pmap, register_openmp
from .analysis.c_action import c_action
from .analysis.c_action import do_action
from .analysis.c_analysis import c_analysis
from .analysis.c_action.actionlist import ActionList
from .topology.topology import Topology
from .builder.build import make_structure
from .analysis.rmsd import (
    rotation_matrix,
    pairwise_rmsd,
    rmsd_perres,
    rmsd_nofit,
    rmsd,
    symmrmsd,
    distance_rmsd, )
from .analysis.energy_analysis import (
    esander,
    lie, )
from .analysis import (
    matrix,
    vector,
    nmr,
    dssp_analysis,
    hbond_analysis,
    energy_analysis, )

from .core.c_core import CpptrajState, Command

# Import all functions from organized modules
from .all_analysis.base_classes import *
from .all_analysis.distance_analysis import *
from .all_analysis.structural_analysis import *
from .all_analysis.geometry_analysis import *
from .all_analysis.transformation import *
from .all_analysis.surface_analysis import *
from .all_analysis.correlation_analysis import *
from .all_analysis.energy_analysis import *
from .all_analysis.struct_utils import *
from .all_analysis.statistical_analysis import *
from .all_analysis.solvation_analysis import *
from .all_analysis.conformational_analysis import *
from .all_analysis.computational_geometry import *

__all__ = [
    'acorr', 'align', 'align_principal_axis', 'analyze_modes', 'angle',
    'atomiccorr', 'atomicfluct', 'atom_map', 'autoimage', 'bfactors',
    'center', 'center_of_geometry', 'center_of_mass', 'check_chirality',
    'check_structure', 'closest', 'closest_atom', 'count_in_voxel', 'crank',
    'density', 'diffusion', 'distance', 'distance_rmsd', 'distance_to_point',
    'distance_to_reference', 'dihedral', 'dssp_analysis', 'energy_analysis',
    'esander', 'fiximagedbonds', 'get_average_frame', 'get_velocity', 'gist',
    'hausdorff', 'hbond_analysis', 'image', 'lie', 'lipidscd', 'lowestcurve',
    'make_structure', 'matrix', 'mean_structure', 'mindist', 'molsurf',
    'multidihedral', 'native_contacts', 'nmr', 'pairdist', 'pairwise_distance',
    'pairwise_rmsd', 'pca', 'permute_dihedrals', 'principal_axes', 'projection',
    'pucker', 'radgyr', 'radgyr_tensor', 'randomize_ions', 'rdf', 'replicate_cell',
    'rmsd', 'rmsd_nofit', 'rmsd_perres', 'rmsf', 'rotate', 'rotate_dihedral',
    'rotation_matrix', 'rotdif', 'scale', 'search_neighbors', 'set_dihedral',
    'set_velocity', 'strip', 'superpose', 'surf', 'symmrmsd', 'ti', 'timecorr',
    'transform', 'translate', 'velocityautocorr', 'vector', 'volmap', 'volume',
    'watershell', 'wavelet', 'xcorr', 'xtalsymm',
] # yapf: disable