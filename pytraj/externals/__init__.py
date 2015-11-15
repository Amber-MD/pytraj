"""this `externals` packages is intended to have
* outside packages
* loading/converting external package's objects to pytraj's object
* loading file format that cpptraj has not yet supported (.h5, ...)
"""
from __future__ import absolute_import
# don't import _load_HDF5 here to avoid circular importing
from ._pickle import to_pickle, read_pickle
from ._json import to_json, read_json

from ._load_MDAnalysis import load_MDAnalysis
from ._load_mdtraj import load_mdtraj
from ._load_ParmEd import load_ParmEd

__all__ = ['read_pickle',
           'read_json',
           'to_pickle',
           'to_json',
           'load_MDAnalysis',
           'load_mdtraj',
           'load_ParmEd', ]
