"""
Notes : might move to cython
"""
from __future__ import absolute_import
from ..datasetlist import DatasetList
from .c_datasetlist import DatasetList as CpptrajDatasetList


def load_datafile(filename):
    """load cpptraj's output"""
    ds = CpptrajDatasetList()
    ds.read_data(filename)
    return DatasetList(ds)
