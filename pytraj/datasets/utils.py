"""
Notes : might move to cython
"""
from __future__ import absolute_import
from ..datasetlist import DatasetList
from ..compat import zip


def load_datafile(filename):
    """load cpptraj's output"""
    ds = DatasetList()
    ds.read_data(filename)
    return ds
