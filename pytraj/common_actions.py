"""having common actions such as rmsd, fitting, ...
>>> from pytraj.common_actions import calc_rmsd, do_fitting
>>> from pytraj.common_actions import translate
# TODO : use __all__
"""
from functools import partial
from pytraj import adict
from pytraj.Topology import Topology
from pytraj.DataSetList import DataSetList
from pytraj.misc import calculate, calc_dssp
from pytraj.DistRoutines import distance 
from pytraj.gdt.calc_score import calc_score

calc_distance = partial(calculate, 'distance')
calc_dih = partial(calculate, 'dihedral')
calc_dihedral = calc_dih
calc_radgyr = partial(calculate, 'radgyr')
calc_angle = partial(calculate, 'angle')
calc_molsurf = partial(calculate, 'molsurf')
calc_molsurf = partial(calculate, 'molsurf')
calc_distrmsd = partial(calculate, 'distrmsd')
calc_protein_score = calc_score

def calc_watershell(command, traj, top=Topology()):
    """return a DataSetList object having the number of water 
    in 1st and 2nd water shell for each frame
    >>> d0 = calc_watershell(":WAT", traj)
    >>> # get 1st shell
    >>> d0_0 = d0[0]
    >>> print (d0_0[:])
    >>> # get 2nd shell
    >>> d0_1 = d0[1]
    >>> print (d0_1[:])
    """
    if 'out' not in command:
        # current Watershell action require specifying output
        # 
        command += ' out .tmp'
    dslist = DataSetList()
    adict['watershell'](command, traj, top, dslist=dslist)
    return dslist
