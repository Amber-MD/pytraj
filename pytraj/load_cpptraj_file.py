"""
for compatibility with cpptraj
>>> from pytraj import load_cpptraj_file
>>> state = load_cpptraj_file(trajin_file)
>>> isinstance(state, CpptrajState)

>>> # state object hold all information like TopologyList, DataSetList,
TrajinList, ...
"""
from pytraj.Command import Command

def load_cpptraj_file(fname):
    """
    Parameters
    ----------
    fname : str, name of cpptraj input file
        ("cpptraj -i input.txt" --> fname = "input.txt")
    """
    return Command.get_state(fname)
