"""
Compatibility module for rmsd functions.
This module re-exports rmsd functions from pytraj.analysis.rmsd
for backward compatibility with code expecting pytraj.actions.rmsd.
"""

# Re-export all rmsd functions from the analysis module
from ..analysis.rmsd import (
    rotation_matrix,
    pairwise_rmsd,
    rmsd_perres,
    rmsd_nofit,
    rmsd,
    symmrmsd,
    distance_rmsd
)

__all__ = [
    'rotation_matrix',
    'pairwise_rmsd',
    'rmsd_perres',
    'rmsd_nofit',
    'rmsd',
    'symmrmsd',
    'distance_rmsd'
]