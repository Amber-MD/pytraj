"""
Notes : might move to cython
"""
from __future__ import absolute_import
from ..datasetlist import DatasetList
from .DataSetList import DataSetList as CpptrajDSL
from ..compat import zip


def load_datafile(filename):
    """load cpptraj's output"""
    ds = CpptrajDSL()
    ds.read_data(filename)
    return DatasetList(ds)
