"""
Compatibility module for matrix functions.
This module re-exports matrix functions from pytraj.analysis.matrix
for backward compatibility with code expecting pytraj.actions.matrix.
"""

# Re-export all matrix functions from the analysis module
from ..analysis.matrix import *