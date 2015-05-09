"""this `externals` packages is intended to have
* outside packages
* loading/converting external package's objects to pytraj's object
* loading file format that cpptraj has not yet supported (.h5, ...)
"""
from __future__ import absolute_import

from ._load_MDAnalysis import  load_MDAnalysis
from ._load_mdtraj import load_mdtraj
from ._load_ParmEd import load_ParmEd
from ._load_pseudo_parm import load_pseudo_parm
