# Import all actions from the new modular structure for backward compatibility

# Import from analysis modules that are still in the original locations
from .analysis.rmsd import (
    rotation_matrix,
    pairwise_rmsd,
    rmsd_perres,
    rmsd_nofit,
    rmsd,
    symmrmsd,
    distance_rmsd
)

from .analysis.energy_analysis import (
    esander,
    lie
)

from .analysis import (
    matrix,
    vector,
    nmr,
    dssp_analysis,
    hbond_analysis,
    energy_analysis
)

from .builder.build import make_structure
from .actions.geometry import _distance_to_ref_or_point
from .actions.correlation import velocity_autocorrelation
from .actions.topology_manipulation import scale as do_scaling
from functools import partial

# Import all actions from the new modular actions package
from .actions import *

# Define the __all__ list for backward compatibility
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
    'watershell', 'wavelet', 'xcorr', 'xtalsymm', 'toroidal_diffusion', 'tordiff',
    'multipucker', 'dihedral_rms', 'ene_decomp', 'infraredspec',
]

# Legacy aliases for backward compatibility
atomicfluct = rmsf
scale = do_scaling
velocityautocorr = velocity_autocorrelation

# Partial functions for distance calculations
distance_to_point = partial(_distance_to_ref_or_point, ref=None)
distance_to_point.__doc__ = """
Compute distance from atoms in mask to a specified point.
Example: pytraj.distance_to_point(traj, ':1', point=[0., 0., 0.])
"""

distance_to_reference = partial(_distance_to_ref_or_point, point=None)
distance_to_reference.__doc__ = """
Compute distance from atoms in mask to a reference structure.
Example: pytraj.distance_to_reference(traj, ':1', ref=ref_frame)
"""

# Helper function for getting average frame
def get_average_frame(traj, mask='*', top=None, **kwargs):
    """Get average structure as a single frame.

    This is a convenience function that returns mean_structure as a single Frame.

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*'
        atom mask
    top : Topology, optional
    **kwargs : additional keyword arguments
        passed to mean_structure (e.g., autoimage, frame_indices)

    Returns
    -------
    Frame
        average frame
    """
    return mean_structure(traj, mask=mask, top=top, **kwargs)

# Alias for toroidal diffusion (same as tordiff)
toroidal_diffusion = tordiff

# Import utility functions
from .actions.base import _assert_mutable, _ensure_mutable
