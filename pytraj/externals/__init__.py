"""this `externals` packages is intended to have
* outside packages
* loading/converting external package's objects to pytraj's object
* loading file format that cpptraj has not yet supported (.h5, ...)
"""
from __future__ import absolute_import
# don't import _load_HDF5 here to avoid circular importing

from .load_other_packages import load_MDAnalysis, load_mdtraj, load_ParmEd

__all__ = [
           'load_MDAnalysis',
           'load_mdtraj',
           'load_ParmEd', ]
