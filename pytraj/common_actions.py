"""having common actions such as rmsd, fitting, ...
>>> from pytraj.common_actions import calc_rmsd
>>> from pytraj.common_actions import translate
# TODO : use __all__
"""
from __future__ import absolute_import
from functools import partial

from pytraj.action_dict import ActionDict
adict = ActionDict()

from pytraj.analysis_dict import AnalysisDict
analdict = AnalysisDict()

from ._common_actions import calculate, _get_top
from .externals.six import string_types
from .Frame import Frame
from .FrameArray import FrameArray
from .AtomMask import AtomMask
from .Topology import Topology
from .DataSetList import DataSetList
from .DataFileList import DataFileList
from .DistRoutines import distance 
from .gdt.calc_score import calc_score
from .hbonds import search_hbonds
from ._shared_methods import _frame_iter_master
from .get_pysander_energies import get_pysander_energies

list_of_cal = ['calc_distance', 'calc_dih', 'calc_dihedral', 'calc_radgyr', 'calc_angle',
               'calc_molsurf', 'calc_distrmsd', 'calc_volume', 'calc_protein_score', 
               'calc_dssp', 'calc_matrix', 'calc_jcoupling',
               'calc_radial', 'calc_watershell',
               'calc_vector',
               'calc_volmap',
               'calc_atomicfluct',
               'calc_COM',
               'calc_center_of_mass',
               'calc_center_of_geometry',
               'calc_pairwise_rmsd']

list_of_do = ['do_translation', 'do_rotation', 'do_autoimage',
              'do_clustering',]

list_of_get = ['get_average_frame']

list_of_the_rest = ['search_hbonds',]

__all__ = list_of_do + list_of_cal + list_of_get + list_of_the_rest

calc_distance = partial(calculate, 'distance', quick_get=True)
calc_dih = partial(calculate, 'dihedral', quick_get=True)
calc_dihedral = calc_dih
calc_radgyr = partial(calculate, 'radgyr', quick_get=True)
calc_angle = partial(calculate, 'angle', quick_get=True)
calc_molsurf = partial(calculate, 'molsurf', quick_get=True)
calc_distrmsd = partial(calculate, 'distrmsd', quick_get=True)
calc_volume = partial(calculate, 'volume', quick_get=True)
calc_matrix = partial(calculate, 'matrix')
calc_jcoupling = partial(calculate, 'jcoupling', quick_get=True)
calc_volmap = partial(calculate, 'volmap', quick_get=True)
calc_protein_score = calc_score
calc_energies = get_pysander_energies

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
    return list(map(lambda idx: ssdict[idx], arr0))

def calc_dssp(command="", traj=None, dtype='int'):
    """return dssp profile for frame/traj

    Parameters
    ---------
    command : str
    traj : {Trajectory, Frame, mix of them}
    dtype : str {'int', 'integer', 'str', 'string', 'dataset'}

    Returns:
    if dtype in ['int', 'integer', 'str', 'string']
        List of tuples with shape (n_frames, n_residues)
    if dtype in ['dataset',]
        DataSetList object
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
    arr0 = list(zip(*arr0))
    if dtype in ['INT', 'INTERGER']:
        return arr0
    elif dtype in ['STRING', 'STR']:
        tmplist = [to_string_ss(arr) for arr in arr0]
        return tmplist
    elif dtype in ['DATASET',]:
        return dslist
    else:
        raise ValueError("")

def do_translation(command="", traj=None, top=Topology()):
    adict['translate'](command, traj, top)

def do_rotation(command="", traj=None, top=Topology()):
    adict['rotate'](command, traj, top)

def do_autoimage(command="", traj=None, top=Topology()):
    adict['autoimage'](command, traj, top)

def get_average_frame(command="", traj=None, top=Topology()):
    dslist = DataSetList()

    # add "crdset s1" to trick cpptraj dumpt coords to DatSetList
    command += " crdset s1"

    act = adict['average']
    act(command, traj, top, dslist=dslist)

    # need to call this method so cpptraj will write
    act.print_output()
    
    return dslist[0].get_frame()

def randomize_ions(command="", traj=Frame(), top=Topology()):
    """randomize_ions for given Frame with Topology
    Return : None
    Parameters
    ---------
    traj : Frame instance, default=Frame()
        frame coords will be modified

    top : Topology instance, default=Topology()

    """
    act = adict['randomizeions']
    act.master(command, traj, top)

def do_clustering(command="", traj=None, top=Topology(), 
        dslist=DataSetList(), dflist=DataFileList()):
    # TODO: still very naive
    ana = analdict['clustering']
    if traj is not None:
        dslist.add_set("coords", "__pytraj_cluster", "")
        dslist[-1].top = top if top is not top.is_empty() else traj.top
        for frame in traj:
            dslist[-1].add_frame(frame)
        command += " crdset __pytraj_cluster"
    if not top.is_empty():
        _top = top
    else:
        _top = traj.top
    ana(command, _top, dslist, dflist) 

def calc_multidihedral(command="", *args, **kwd): 
    """perform dihedral search
    Parameters
    ----------
    command : str, cpptraj command 
    traj : Trajectory-like object
    *arg and **kwd: additional arguments

    Returns
    -------
    Dictionary of array
    >>> from pytraj.common_actions import calc_multidihedral
    >>> d = calc_multidihedral("resrange 6-9 phi psi chi", traj)
    >>> assert isinstance(d, dict) == True
    >>> from pytraj.dataframe import to_dataframe
    >>> print (to_dataframe(d))
    """
    dslist = DataSetList()
    from pytraj.six_2 import izip as zip
    from array import array
    act = adict['multidihedral']
    act(command, dslist=dslist, *args, **kwd)
    return dict((d0.legend, array('d', d0.data)) for d0 in dslist)

def calc_atomicfluct(command="", *args, **kwd):
    dslist = DataSetList()
    dslist.set_py_free_mem(False)
    act = adict['atomicfluct']
    act(command, dslist=dslist, *args, **kwd)
    # tag: print_output()
    act.print_output() # need to have this. check cpptraj's code
    return dslist[-1]

def calc_vector(mask="", traj=None, *args, **kwd): 
    """perform dihedral search
    Parameters
    ----------
    command : str, cpptraj command 
    traj : Trajectory-like object
    *arg and **kwd: additional arguments

    Returns
    -------
    DataSet_Vector object

    Examples
    ------
    >>> import pytraj.common_actions as pyca
    >>> pyca.calc_vector("@CA @CB", traj).tolist()
    >>> pyca.calc_vector("", traj).tolist()
    >>> pyca.calc_vector("principal z", traj).to_ndarray()
    >>> pyca.calc_vector("principal x", traj).to_ndarray()
    >>> pyca.calc_vector("ucellx", traj).tolist()
    >>> pyca.calc_vector("boxcenter", traj).tolist()
    >>> pyca.calc_vector("box", traj).tolist()
    """
    from pytraj.actions.Action_Vector import Action_Vector
    from pytraj.DataSetList import DataSetList
    act = Action_Vector()
    dslist = DataSetList()

    if 'name' not in mask:
        # for some reasons, I got segmentation fault without 'name' keyword
        # need to check cpptraj code
        mask = "name myvector " + mask
    act(command=mask, current_frame=traj, dslist=dslist, *args, **kwd)
    dslist.set_py_free_mem(False)
    return dslist[0]

def _calc_vector_center(command="", traj=None, top=None, use_mass=False):
    _top = _get_top(traj, top)

    dslist = DataSetList()
    dslist.set_py_free_mem(False) # need this to avoid segmentation fault
    act = adict['vector']
    command = "center " + command
    act.read_input(command=command, current_top=_top, dslist=dslist)
    act.process(_top)
    for frame in _frame_iter_master(traj):
        # set Frame masses
        if use_mass:
            frame.set_frame_m(_top)
        act.do_action(frame)
    return dslist[0]

calc_COM = calc_center_of_mass = partial(_calc_vector_center, use_mass=True)
calc_COG = calc_center_of_geometry = partial(_calc_vector_center, use_mass=False)

def calc_pairwise_rmsd(command="", traj=None, top=None, *args, **kwd):
    """return  DataSetList object

    Examples:
        dslist = calc_pairwise_rmsd("@CA", traj)
        dslist.to_ndarray()
        dslist.tolist()
    """
    from pytraj.analyses import Analysis_Rms2d
    act = Analysis_Rms2d()

    dslist = DataSetList()
    dslist.add_set("coords", "mylovelypairwisermsd")
    _top = _get_top(traj, top)
    dslist[0].top = _top
    # need to set "rmsout" to trick cpptraj not giving error
    # need " " (space) before crdset too
    command = command + " crdset mylovelypairwisermsd rmsout mycrazyoutput"

    # upload Frame to crdset
    for frame in _frame_iter_master(traj):
        dslist[0].append(frame)

    act(command, _top, dslist=dslist, *args, **kwd)
    # remove dataset coords to free memory
    dslist.remove_set(dslist[0])
    return dslist
