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
