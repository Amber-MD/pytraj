"""
for compatibility with cpptraj
>>> from pytraj import load_cpptraj_file
>>> state = load_cpptraj_file(trajin_file)
>>> isinstance(state, CpptrajState)

>>> # state object hold all information like TopologyList, DataSetList,
TrajinList, ...
"""
from __future__ import absolute_import
from .core.Command import Command
from .utils import file_exist as file_exists

def load_cpptraj_file(fname):
    """
    Parameters
    ----------
    fname : str, name of cpptraj input file
        ("cpptraj -i input.txt" --> fname = "input.txt")
    """
    if not file_exists(fname):
        raise ValueError("can not locate this file")
    return Command.get_state(fname)
