"""
for compatibility with cpptraj
>>> from pytraj import load_cpptraj_file
>>> state = load_cpptraj_file(trajin_file)
>>> isinstance(state, CpptrajState)
"""
from __future__ import absolute_import
from ..core.c_core import Command
from ..decorators import ensure_exist


@ensure_exist
def load_cpptraj_file(filename):
    """
    Parameters
    ----------
    fname : str, name of cpptraj input file
        ("cpptraj -i input.txt" --> fname = "input.txt")
    """
    with Command() as cm:
        return cm.get_state(filename)
