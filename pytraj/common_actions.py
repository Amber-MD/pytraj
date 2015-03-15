"""having common actions such as rmsd, fitting, ...
>>> from pytraj.common_actions import calc_rmsd
>>> from pytraj.common_actions import translate
# TODO : use __all__
"""
from pytraj._common_actions import calculate
from functools import partial
from pytraj import adict
from pytraj.externals.six import string_types
from pytraj.Frame import Frame
from pytraj.FrameArray import FrameArray
from pytraj.AtomMask import AtomMask
from pytraj.Topology import Topology
from pytraj.DataSetList import DataSetList
from pytraj.DistRoutines import distance 
from pytraj.gdt.calc_score import calc_score

list_of_cal = ['calc_distance', 'calc_dih', 'calc_dihedral', 'calc_radgyr', 'calc_angle',
               'calc_molsurf', 'calc_distrmsd', 'calc_volume', 'calc_protein_score', 
               'calc_dssp', 'calc_dssp',
               'calc_watershell']

list_of_do = ['do_translation', 'do_rotation', 'do_autoimage']

__all__ = list_of_do + list_of_cal

calc_distance = partial(calculate, 'distance')
calc_dih = partial(calculate, 'dihedral')
calc_dihedral = calc_dih
calc_radgyr = partial(calculate, 'radgyr')
calc_angle = partial(calculate, 'angle')
calc_molsurf = partial(calculate, 'molsurf')
calc_distrmsd = partial(calculate, 'distrmsd')
calc_volume = partial(calculate, 'volume')
calc_protein_score = calc_score

do_translation = partial(calculate, 'translate')
do_rotation = partial(calculate, 'rotate')

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

def calc_radial(command, traj, top=Topology()):
    '''Action_Radial require calling Print() to get output. We make change here'''
    act = adict['radial']
    # add `radial` keyword to command (need to check `why`?)
    command = 'radial ' + command
    dslist = DataSetList()
    if not top.is_empty():
        act(command, traj, top, dslist=dslist)
    else:
        act(command, traj, dslist=dslist)

    # dump data to dslist.
    act.print_output()
    return dslist

def to_string_ss(arr0):
    """
    arr0 : ndarray
    """
    #ss = ['None', 'Para', 'Anti', '3-10', 'Alpha', 'Pi', 'Turn', 'Bend']
    ss = ["0", "b", "B", "G", "H", "I", "T", "S"]
    len_ss = len(ss)
    ssdict = dict(zip(range(len_ss), ss))
    return map(lambda idx: ssdict[idx], arr0)

def calc_dssp(command="", traj=None, dtype='int'):
    """return dssp profile for frame/traj

    Parameters
    ---------
    command : str
    traj : {Trajectory, Frame, mix of them}
    dtype : str {'int', 'integer', 'str', 'string'}

    Returns:
    List of tuples with shape (n_frames, n_residues)
    """
    dslist = DataSetList()
    adict['dssp'](command,
                  current_frame=traj, current_top=traj.top,
                  dslist=dslist)
    dtype = dtype.upper()

    # get all dataset from DatSetList if dtype == integer
    arr0 = list(dslist.get_dataset(dtype="integer"))

    # cpptraj store data for each residue for each frame (n_residues, n_frames)
    # we need to transpose data
    arr0 = zip(*arr0)
    if dtype in ['INT', 'INTERGER']:
        return arr0
    elif dtype in ['STRING', 'STR']:
        tmplist = [to_string_ss(arr) for arr in arr0]
        return tmplist

def do_translation(command="", traj=None, top=Topology()):
    adict['translate'](command, traj, top)

def do_rotation(command="", traj=None, top=Topology()):
    adict['rotate'](command, traj, top)

def do_autoimage(command="", traj=None, top=Topology()):
    adict['autoimage'](command, traj, top)

