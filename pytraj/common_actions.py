"""having common actions such as rmsd, fitting, ...
>>> from pytraj.common_actions import calc_rmsd
>>> from pytraj.common_actions import translate
# TODO : use __all__
"""
from __future__ import absolute_import
from array import array
from functools import partial

from pytraj.action_dict import ActionDict
adict = ActionDict()

from pytraj.analysis_dict import AnalysisDict
analdict = AnalysisDict()

from ._get_common_objects import _get_top, _get_data_from_dtype
from ._common_actions import calculate
from .externals.six import string_types
from .Frame import Frame
from .FrameArray import FrameArray
from .AtomMask import AtomMask
from .Topology import Topology
from .DataSetList import DataSetList
from .DataFileList import DataFileList
from .math.DistRoutines import distance 
from .externals.gdt.calc_score import calc_score
from .hbonds import search_hbonds
from ._shared_methods import _frame_iter_master
from .externals.get_pysander_energies import get_pysander_energies
from .utils import _import_numpy, is_array

list_of_cal = ['calc_distance', 'calc_dih', 'calc_dihedral', 'calc_radgyr', 'calc_angle',
               'calc_molsurf', 'calc_distrmsd', 'calc_volume', 'calc_protein_score', 
               'calc_dssp', 'calc_matrix', 'calc_jcoupling',
               'calc_radial', 'calc_watershell',
               'calc_vector',
               'calc_multivector',
               'calc_volmap',
               'calc_rdf',
               'calc_atomicfluct',
               'calc_COM',
               'calc_center_of_mass',
               'calc_center_of_geometry',
               'calc_pairwise_rmsd',
               'calc_temperatures']

list_of_do = ['do_translation', 'do_rotation', 'do_autoimage',
              'do_clustering',]

list_of_get = ['get_average_frame']

list_of_the_rest = ['search_hbonds', 'align_principal_axis', 'closest']

__all__ = list_of_do + list_of_cal + list_of_get + list_of_the_rest

#calc_distance = partial(calculate, 'distance', quick_get=True)
calc_dih = partial(calculate, 'dihedral', quick_get=True)
calc_dihedral = calc_dih
calc_radgyr = partial(calculate, 'radgyr', quick_get=True)
calc_angle = partial(calculate, 'angle', quick_get=True)
calc_molsurf = partial(calculate, 'molsurf', quick_get=True)
calc_distrmsd = partial(calculate, 'distrmsd', quick_get=True)
calc_volume = partial(calculate, 'volume', quick_get=True)
calc_matrix = partial(calculate, 'matrix')
calc_jcoupling = partial(calculate, 'jcoupling', quick_get=True)
calc_multivector = partial(calculate, 'multivector')
calc_volmap = partial(calculate, 'volmap', quick_get=True)
calc_rdf = partial(calculate, 'radial', print_output=True)
calc_protein_score = calc_score
calc_energies = get_pysander_energies

action_type = calculate
do_translation = partial(action_type, 'translate')
translate = do_translation
do_rotation = partial(action_type, 'rotate')
rotate = do_rotation
do_scaling = partial(action_type, 'scale')
scale = do_scaling

def calc_distance(traj=None, command="", top=None, *args, **kwd):
    """calculate distance

    Notes:
    command : str | int_2d numpy array
    """
    _, np = _import_numpy()
    _top = _get_top(traj, top)
    if isinstance(command, string_types):
        # cpptraj mask for action
        return calculate("distance", traj, command, top=_top, quick_get=True, *args, **kwd)
    elif isinstance(command, np.ndarray):
        int_2darr = command
        if 'n_frames' not in kwd.keys():
            raise ValueError("require specifying n_frames")
        arr = np.empty([kwd['n_frames'], len(int_2darr)])
        for idx, frame in enumerate(_frame_iter_master(traj)):
            arr[idx] = frame.calc_distance(int_2darr)
        return arr
    else:
        raise ValueError("")

def calc_mindist(traj=None, command="", top=None, *args, **kwd):
    _command = "mindist " + command 
    _top = _get_top(traj, top)
    return calculate("nativecontacts", traj, _command, top=_top, quick_get=True, *args, **kwd)

def calc_watershell(traj=None, command="", top=Topology()):
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
    _top = _get_top(traj, top)
    if not 'out' in command:
        # current Watershell action require specifying output
        # 
        command += ' out .tmp'
    dslist = DataSetList()
    adict['watershell'](command, traj, _top, dslist=dslist)
    return dslist

def calc_radial(traj=None, command="", top=Topology()):
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

def calc_dssp(traj=None, command="", top=None, dtype='int'):
    """return dssp profile for frame/traj

    Parameters
    ---------
    command : str
    traj : {Trajectory, Frame, mix of them}
    dtype : str {'int', 'integer', 'str', 'string', 'dataset', 'ndarray'}

    Returns:
    if dtype in ['int', 'integer', 'str', 'string']
        List of tuples with shape (n_frames, n_residues)
    if dtype in ['dataset',]
        DataSetList object
    """
    _top = _get_top(traj, top)
    dslist = DataSetList()
    adict['dssp'](command,
                  current_frame=traj, current_top=_top,
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
    if dtype in ['NDARRAY',]:
        # return a numpy array of strings
        _, np = _import_numpy()
        return np.array([to_string_ss(arr) for arr in arr0])
    else:
        raise ValueError("")

def do_translation(traj=None, command="", top=Topology()):
    adict['translate'](command, traj, top)

def do_rotation(traj=None, command="",  top=Topology()):
    adict['rotate'](command, traj, top)

def do_autoimage(traj=None, command="", top=Topology()):
    adict['autoimage'](command, traj, top)

autoimage = do_autoimage

def get_average_frame(traj=None, command="", top=Topology()):
    dslist = DataSetList()

    # add "crdset s1" to trick cpptraj dumpt coords to DatSetList
    command += " crdset s1"

    act = adict['average']
    act(command, traj, top, dslist=dslist)

    # need to call this method so cpptraj will write
    act.print_output()
    
    return dslist[0].get_frame()

def randomize_ions(traj=Frame(), command="", top=Topology()):
    """randomize_ions for given Frame with Topology
    Return : None
    Parameters
    ---------
    traj : Frame instance, default=Frame()
        frame coords will be modified

    top : Topology instance, default=Topology()

    """
    act = adict['randomizeions']
    act(command, traj, top)

def do_clustering(traj=None, command="", top=Topology(), 
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

def calc_multidihedral(traj=None, command="", dtype='dict', top=None, *args, **kwd): 
    """perform dihedral search
    Parameters
    ----------
    command : str, cpptraj command 
    traj : Trajectory-like object
    *arg and **kwd: additional arguments

    Returns
    -------
    Dictionary of array or dataset or ndarray or list or pyarray (based on `dtype`)

    Examples
    --------
        from pytraj.common_actions import calc_multidihedral
        # calculate all phi/psi dihedrals for residues 6 to 9 (for Amber string mask, index starts from 1)
        d = calc_multidihedral(traj, "resrange 6-9 phi psi chi")
        assert isinstance(d, dict) == True
        from pytraj.dataframe import to_dataframe
        print (to_dataframe(d))

        # calculate dihedrals for N:CA:CB:CG for all residues, return 'DataSetList' object 
        d = calc_multidihedral(traj, "dihtype chi1:N:CA:CB:CG", dtype='dataset'))

        # calculate all dihedrals, save output to DataSetList and write output to disk too 
        from pytraj import DataFileList
        dflist = DataFileList()
        d = pdb.calc_multidihedral("out ./output/test_multdih.dat", dtype='dataset', dflist=dflist)
        dflist.write_all_datafiles()

    See Also
    -------
        Amber15 manual: http://ambermd.org/doc12/Amber15.pdf (page 579)
    """
    _top = _get_top(traj, top)
    dslist = DataSetList()
    from pytraj.six_2 import izip as zip
    from array import array
    act = adict['multidihedral']
    act(command, traj, _top, dslist=dslist, *args, **kwd)
    if dtype == 'dict':
        return dict((d0.legend, array('d', d0.data)) for d0 in dslist)
    else:
        # return dslist
        return dslist

def calc_atomicfluct(traj=None, command="", *args, **kwd):
    dslist = DataSetList()
    dslist.set_py_free_mem(False)
    act = adict['atomicfluct']
    act(command, traj, dslist=dslist, *args, **kwd)
    # tag: print_output()
    act.print_output() # need to have this. check cpptraj's code
    return dslist[-1]

def calc_vector(traj=None, mask="", dtype='vector', *args, **kwd): 
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
    >>> pyca.calc_vector(traj, "@CA @CB").tolist()
    >>> pyca.calc_vector(traj, "", traj).tolist()
    >>> pyca.calc_vector(traj, "principal z").to_ndarray()
    >>> pyca.calc_vector(traj, "principal x").to_ndarray()
    >>> pyca.calc_vector(traj, "ucellx").tolist()
    >>> pyca.calc_vector(traj, "boxcenter").tolist()
    >>> pyca.calc_vector(traj, "box").tolist()
    """
    from pytraj.actions.Action_Vector import Action_Vector
    from pytraj.DataSetList import DataSetList
    act = Action_Vector()
    dslist = DataSetList()

    act(command=mask, current_frame=traj, dslist=dslist, *args, **kwd)
    dslist.set_py_free_mem(False)
    if dtype == 'vector':
        return dslist[-1]
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)

def _calc_vector_center(traj=None, command="", top=None, use_mass=False, dtype='dataset'):
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
    return _get_data_from_dtype(dslist[0], dtype=dtype)

calc_COM = calc_center_of_mass = partial(_calc_vector_center, use_mass=True)
calc_COG = calc_center_of_geometry = partial(_calc_vector_center, use_mass=False)

def calc_pairwise_rmsd(traj=None, command="", top=None, *args, **kwd):
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

def calc_temperatures(traj=None, command="", top=None):
    """return 1D python array of temperatures (from velocity) in traj
    if `frame` keyword is specified cpptraj/pytraj will take existing T

    Default = array of 0.0
    """
    from array import array as pyarray
    _top = _get_top(traj, top)
    dslist = calculate('temperature', traj, command, _top)
    return pyarray('d', dslist[0].tolist())

def calc_rmsd(traj=None, command="", ref=None, mass=False, fit=True, top=None):
    """calculate rmsd

    Parameters
    ---------
    command : str
        Atom mask
    traj : Trajectory | List of trajectories | Trajectory or frame_iter
    top : Topology | str
        (optional) Topology
    ref : Frame | str, default=None (ust 1st frame)
    mass : bool, default=True
        use mass or not
    fit : bool, default=True
        fit or no fit

    Examples
    --------
        calc_rmsd(traj, ":3-18@CA", ref=traj[0], mass=True, fit=True)
        calc_rmsd(traj, ":3-18@CA", 0) # ref=traj[0]
        calc_rmsd(traj, ":3-18@CA", 'last') # ref=traj[-1]
        calc_rmsd(traj, ":3-18@CA", 'first') # ref=traj[0]
        calc_rmsd(traj, ":3-18@CA", 'Tc5b.nat.crd') # ref: from file

    """
    _top = _get_top(traj, top)
    if ref is None or ref == 'first':
        # set ref to 1st frame
        ref = traj[0]
    elif ref == 'last':
        ref = traj[-1]
    elif isinstance(ref, string_types):
        # need to check this in the end to avoid using 'last' keyword
        from .trajs.Trajin_Single import Trajin_Single
        ref = Trajin_Single(ref, _top)[0]

    arr = array('d')

    # creat AtomMask object
    if isinstance(command, string_types):
        atm = _top(command) 
    elif isinstance(command, AtomMask):
        atm = command
    elif is_array(command) or isinstance(command, (list, tuple)):
        atm = AtomMask()
        atm.add_selected_indices(command)
    else:
        atm = AtomMask()

    if mass:
        ref.set_frame_m(_top)
    _ref = Frame(ref, atm)
    for frame in traj:
        if mass:
            # TODO : just need to set mass once
            frame.set_frame_m(_top)
        if fit:
            _rmsd = frame.rmsd(ref, atommask=atm, use_mass=mass)
        else:
            _frame = Frame(frame, atm)
            _rmsd = _frame.rmsd_nofit(ref, use_mass=mass)
        arr.append(_rmsd)
    return arr

def align_principal_axis(traj=None, command="*", top=None):
    # TODO : does not match with cpptraj output
    # rmsd_nofit ~ 0.5 for md1_prod.Tc5b.x, 1st frame
    """
    Notes
    -----
    apply for mutatble traj (FrameArray, Frame)
    """
    act = adict['principal']
    command += " dorotation"
    act(command, traj, top)

def closest(traj=None, command=None, dslist=None, top=None, *args, **kwd):
    from .actions.Action_Closest import Action_Closest
    act = Action_Closest()
    fa = FrameArray()

    _top = _get_top(traj, top)
    new_top = Topology()
    new_top.py_free_mem = False # cpptraj will do
    act.read_input(command, _top)
    act.process(_top, new_top)

    fa.top = new_top.copy()
    for frame in _frame_iter_master(traj):
        new_frame = Frame()
        new_frame.py_free_mem = False # cpptraj will do
        act.do_action(frame, new_frame)
        fa.append(new_frame, copy=True)
    return fa
