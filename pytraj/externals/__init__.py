"""this `externals` packages is intended to have
* outside packages
* loading/converting external package's objects to pytraj's object
* loading file format that cpptraj has not yet supported (.h5, ...)
"""
from .load_other_packages import load_parmed

__all__ = [
    'load_parmed',
]
