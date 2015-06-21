"""having common actions such as rmsd, fitting, ...
>>> from pytraj.common_actions import calc_rmsd
>>> from pytraj.common_actions import translate
# TODO : use __all__
"""
from __future__ import absolute_import
import os
from array import array
from functools import partial

from pytraj.action_dict import ActionDict
adict = ActionDict()

from pytraj.analysis_dict import AnalysisDict
analdict = AnalysisDict()

from ._get_common_objects import _get_top, _get_data_from_dtype, _get_list_of_commands
from ._common_actions import calculate
from .utils import _import_numpy, is_array, ensure_not_none_or_string
from .externals.six import string_types
from .Frame import Frame
#from .Trajectory import Trajectory
from .AtomMask import AtomMask
from .Topology import Topology
from .datasets.DataSetList import DataSetList
from .core.DataFileList import DataFileList
from .math.DistRoutines import distance 
from .externals.gdt.calc_score import calc_score
from .hbonds import search_hbonds, search_nointramol_hbonds
from ._shared_methods import _frame_iter_master
from .externals.get_pysander_energies import get_pysander_energies
from .utils.context import goto_temp_folder
from . import _long_manual

list_of_cal = ['calc_distance', 'calc_dihedral', 'calc_radgyr', 'calc_angle',
               'calc_molsurf', 'calc_distrmsd', 'calc_volume', 'calc_protein_score', 
               'calc_dssp', 'calc_matrix', 'calc_jcoupling',
               'calc_radial', 'calc_watershell',
               'calc_vector',
               'calc_multivector',
               'calc_volmap',
               'calc_rdf',
               'calc_multidihedral',
               'calc_atomicfluct',
               'calc_COM',
               'calc_center_of_mass',
               'calc_center_of_geometry',
               'calc_pairwise_rmsd',
               'calc_density',
               'calc_grid',
               'calc_temperatures',
               'calc_linear_interaction_energy',]

list_of_do = ['do_translation', 'do_rotation', 'do_autoimage',
              'do_clustering',]

list_of_get = ['get_average_frame']

list_of_the_rest = ['search_hbonds', 'search_nointramol_hbonds', 
                    'align_principal_axis', 'closest',
                    'native_contacts', 'nastruct']

__all__ = list_of_do + list_of_cal + list_of_get + list_of_the_rest

calc_protein_score = calc_score
calc_energies = get_pysander_energies
energy_decomposition = get_pysander_energies

action_type = calculate
do_translation = partial(action_type, 'translate')
translate = do_translation
do_rotation = partial(action_type, 'rotate')
rotate = do_rotation
do_scaling = partial(action_type, 'scale')
scale = do_scaling

def calc_distance(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    """calculate distance

    Notes:
    command : str | list of strings | int_2d numpy array
    """
    ensure_not_none_or_string(traj)

    _, np = _import_numpy()
    _top = _get_top(traj, top)
    if isinstance(command, string_types):
        # need to remove 'n_frames' keyword since Action._master does not use it
        try:
            del kwd['n_frames']
        except:
            pass
        # cpptraj mask for action
        dset = calculate("distance", traj, command, top=_top,  quick_get=True, *args, **kwd)
        return _get_data_from_dtype(dset, dtype)
    elif isinstance(command, (list, tuple)):
        list_of_commands = command
        from pytraj.core.ActionList import ActionList
        from pytraj.actions.CpptrajActions import Action_Distance
        dslist = DataSetList()
        actlist = ActionList()

        for cm in list_of_commands:
            actlist.add_action(Action_Distance(), cm, _top, dslist=dslist, *args, **kwd)
        actlist.do_actions(traj)
        return _get_data_from_dtype(dslist, dtype)
    elif isinstance(command, np.ndarray):
        int_2darr = command
        if int_2darr.shape[1]  != 2:
            raise ValueError("require int-array with shape=(n_atoms, 2)")
        if 'n_frames' not in kwd.keys():
            try:
                n_frames = traj.n_frames
            except:
                raise ValueError("require specifying n_frames")
        else:
            n_frames = kwd['n_frames']
        arr = np.empty([n_frames, len(int_2darr)])
        for idx, frame in enumerate(_frame_iter_master(traj)):
            arr[idx] = frame.calc_distance(int_2darr)
        return arr
    else:
        raise ValueError("")

def calc_angle(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    """calculate dihedral

    Notes:
    command : str | int_2d numpy array
    """
    ensure_not_none_or_string(traj)

    _, np = _import_numpy()
    _top = _get_top(traj, top)
    if isinstance(command, string_types):
        # need to remove 'n_frames' keyword since Action._master does not use it
        try:
            del kwd['n_frames']
        except:
            pass
        # cpptraj mask for action
        dset = calculate("angle", traj, command, top=_top, quick_get=True, *args, **kwd)
        return _get_data_from_dtype(dset, dtype)
    elif isinstance(command, (list, tuple)):
        list_of_commands = command
        from pytraj.core.ActionList import ActionList
        from pytraj.actions.CpptrajActions import Action_Angle
        dslist = DataSetList()
        actlist = ActionList()

        for cm in list_of_commands:
            actlist.add_action(Action_Angle(), cm, _top, dslist=dslist, *args, **kwd)
        actlist.do_actions(traj)
        return _get_data_from_dtype(dslist, dtype)
    elif isinstance(command, np.ndarray):
        int_2darr = command
        if int_2darr.shape[1]  != 3:
            raise ValueError("require int-array with shape=(n_atoms, 3)")
        if 'n_frames' not in kwd.keys():
            try:
                n_frames = traj.n_frames
            except:
                raise ValueError("require specifying n_frames")
        else:
            n_frames = kwd['n_frames']
        arr = np.empty([n_frames, len(int_2darr)])
        for idx, frame in enumerate(_frame_iter_master(traj)):
            arr[idx] = frame.calc_angle(int_2darr)
        return arr
    else:
        raise ValueError("")

def calc_dihedral(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    """calculate dihedral

    Notes:
    command : str | int_2d numpy array
    """
    ensure_not_none_or_string(traj)

    _, np = _import_numpy()
    _top = _get_top(traj, top)
    if isinstance(command, string_types):
        # need to remove 'n_frames' keyword since Action._master does not use it
        try:
            del kwd['n_frames']
        except:
            pass
        # cpptraj mask for action
        dset = calculate("dihedral", traj, command, top=_top, quick_get=True, *args, **kwd)
        return _get_data_from_dtype(dset, dtype)
    elif isinstance(command, (list, tuple)):
        list_of_commands = command
        from pytraj.core.ActionList import ActionList
        from pytraj.actions.CpptrajActions import Action_Dihedral
        dslist = DataSetList()
        actlist = ActionList()

        for cm in list_of_commands:
            actlist.add_action(Action_Dihedral(), cm, _top, dslist=dslist, *args, **kwd)
        actlist.do_actions(traj)
        return _get_data_from_dtype(dslist, dtype)
    elif isinstance(command, np.ndarray):
        int_2darr = command
        if int_2darr.shape[1]  != 4:
            raise ValueError("require int-array with shape=(n_atoms, 4)")
        if 'n_frames' not in kwd.keys():
            try:
                n_frames = traj.n_frames
            except:
                raise ValueError("require specifying n_frames")
        else:
            n_frames = kwd['n_frames']
        arr = np.empty([n_frames, len(int_2darr)])
        for idx, frame in enumerate(_frame_iter_master(traj)):
            arr[idx] = frame.calc_dihedral(int_2darr)
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

def calc_matrix(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_Matrix
    act = Action_Matrix()

    _top = _get_top(traj, top)
    dslist = DataSetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    act.print_output()
    return _get_data_from_dtype(dslist, dtype)


def calc_radgyr(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_Radgyr
    act = Action_Radgyr()

    _top = _get_top(traj, top)
    dslist = DataSetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist[0].copy(), dtype)


def calc_molsurf(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_Molsurf
    act = Action_Molsurf()

    _top = _get_top(traj, top)
    dslist = DataSetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_distrmsd(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_DistRmsd
    act = Action_DistRmsd()

    _top = _get_top(traj, top)
    dslist = DataSetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_volume(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_Volume
    act = Action_Volume()

    _top = _get_top(traj, top)
    dslist = DataSetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_multivector(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_MultiVector
    act = Action_MultiVector()

    _top = _get_top(traj, top)
    dslist = DataSetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_volmap(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_Volmap
    act = Action_Volmap()

    _top = _get_top(traj, top)
    dslist = DataSetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_linear_interaction_energy(traj=None, command="", top=None, 
        dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_LIE
    act = Action_LIE()

    _top = _get_top(traj, top)
    dslist = DataSetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)

# alias
calc_LIE = calc_linear_interaction_energy


def calc_rdf(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_Radial
    act = Action_Radial()

    _top = _get_top(traj, top)
    dslist = DataSetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    act.print_output()
    return _get_data_from_dtype(dslist, dtype)


def calc_jcoupling(traj=None, command="", top=None, kfile=None, dtype='dataset', *args, **kwd):
    """
    Paramters
    ---------
    traj : any things that make `frame_iter_master` returning Frame object
    command : str, default ""
        cpptraj's command/mask
    kfile : str, default None, optional
        Dir for Karplus file. If "None", use $AMBERHOME dir 
    dtype : str, {'dataset', ...}, default 'dataset'
    *args, **kwd: optional
    """
    from pytraj.actions.CpptrajActions import Action_Jcoupling
    act = Action_Jcoupling()
    # add `radial` keyword to command (need to check `why`?)
    dslist = DataSetList()
    _top = _get_top(traj, top)

    if kfile is not None:
        command += " kfile %s" % kfile
    act(command, traj, dslist=dslist, top=_top, *args, **kwd)
    return dslist

def to_string_ss(arr0):
    """
    arr0 : ndarray
    """
    _, np = _import_numpy()
    #ss = ['None', 'Para', 'Anti', '3-10', 'Alpha', 'Pi', 'Turn', 'Bend']
    ss = ["0", "b", "B", "G", "H", "I", "T", "S"]
    len_ss = len(ss)
    ssdict = dict(zip(range(len_ss), ss))

    if np:
        def myfunc(key):
            return ssdict[key]
        if not isinstance(arr0, dict):
            return np.vectorize(myfunc)(arr0)
        else:
            new_dict = {}
            for key in arr0.keys():
                new_dict[key] = to_string_ss(arr0[key])
            return new_dict
    else:
        print ("doest not have numpy, return a list")
        return list(map(lambda idx: ssdict[idx], arr0))

def calc_dssp(traj=None, command="", top=None, dtype='int', dslist=None, dflist=DataFileList()):
    """return dssp profile for frame/traj

    Parameters
    ----------
    command : str
    traj : {Trajectory, Frame, mix of them}
    dtype : str {'int', 'integer', 'str', 'string', 'dataset', 'ndarray'}

    Returns
    -------
    if dtype in ['int', 'integer', 'str', 'string']
        List of tuples with shape (n_frames, n_residues)
    if dtype in ['dataset',]
        DataSetList object

    Examples
    --------
        calc_dssp(traj, ":2-10")

        calc_dssp(traj, ":2-10 out dssp.gnu", dflist=dflist)
        dflist.write_all_datafiles()

        calc_dssp(traj, ":2-10 sumout dssp.agr", dflist=dflist)
        dflist.write_all_datafiles()
        # from terminal: xmgrace dssp.agr

    Notes
    -----
    Character Integer DSSP_Char SS_type
    0         0       ' '       None
    b         1       'E'       Parallel Beta-sheet
    B         2       'B'       Anti-parallel Beta-sheet
    G         3       'G'       3-10 helix
    H         4       'H'       Alpha helix
    I         5       'I'       Pi (3-14) helix
    T         6       'T'       Turn
    S         7       'S'       Bend

    See Also
    --------
    Amber15 manual: http://ambermd.org/doc12/Amber15.pdf (page 588)
    """
    _, np = _import_numpy()

    _top = _get_top(traj, top)
    if dslist is None:
        dslist = DataSetList()
    adict['dssp'](command,
                  current_frame=traj, 
                  top=_top,
                  dslist=dslist,
                  dflist=dflist)

    # replace legend to something nicer
    for legend, dset in dslist.iteritems():
        if 'DSSP' in legend:
            legend = legend.replace("DSSP_00000[", "")
            legend = legend.replace("]", "_avg")
            dset.legend = legend.lower()
    dtype = dtype.upper()

    # get all dataset from DatSetList if dtype == integer
    if not np:
        arr0 = list(dslist.get_dataset(dtype="integer"))

        # cpptraj store data for each residue for each frame (n_residues, n_frames)
        # we need to transpose data
        arr0 = list(zip(*arr0))
    else:
        arr0 = dslist.groupby("integer", mode='dtype').values
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
        try:
            return _get_data_from_dtype(dslist, dtype)
        except:
            raise ValueError("")

def do_translation(traj=None, command="", top=Topology()):
    adict['translate'](command, traj, top)

def do_rotation(traj=None, command="",  top=Topology()):
    adict['rotate'](command, traj, top)

def do_autoimage(traj=None, command="", top=Topology()):
    adict['autoimage'](command, traj, top)

autoimage = do_autoimage

def get_average_frame(traj=None, command="", top=Topology()):
    _top = _get_top(traj, top)
    dslist = DataSetList()

    # add "crdset s1" to trick cpptraj dumpt coords to DatSetList
    command += " crdset s1"

    act = adict['average']
    act(command, traj, _top, dslist=dslist)

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

def do_clustering(traj=None, command="", top=None, dtype='dataset',
        dslist=None, dflist=None):
    """
    Parameters
    ---------
    traj : Trajectory-like | list of Trajectory-like | frame or chunk iterator
    command : cpptraj command
    top : Topology, optional
    dslist : DataSetList, optional
    dflist : DataFileList, optional

    Notes:
    Supported algorithms: kmeans, hieragglo, and dbscan.

    Examples
    --------
        do_clustering(traj, "kmeans clusters 50 @CA")

    Returns
    -------
    DataSetList object
 
    """

    _top = _get_top(traj, top)
    ana = analdict['clustering']
    # need to creat `dslist` here so that every time `do_clustering` is called,
    # we will get a fresh one (or will get segfault)
    if dslist is None:
        dslist = DataSetList()
    else:
        dslist = dslist

    if traj is not None:
        dslist.add_set("coords", "__pytraj_cluster")
        #dslist[-1].top = _top
        dslist[0].top = _top
        for frame in traj:
            #dslist[-1].add_frame(frame)
            dslist[0].add_frame(frame)
        command += " crdset __pytraj_cluster"
    else:
        pass
    ana(command, _top, dslist, dflist) 
    # remove frames in dslist to save memory
    dslist.remove_set(dslist['__pytraj_cluster'])
    return _get_data_from_dtype(dslist, dtype=dtype)

def calc_multidihedral(traj=None, command="", dtype='dataset', top=None, *args, **kwd): 
    """perform dihedral search
    Parameters
    ----------
    command : str, cpptraj command 
    traj : Trajectory-like object
    *arg and **kwd: additional arguments

    Returns
    -------
    Dictionary of array or dataset or ndarray or list or pyarray (based on `dtype`)

    Notes
    -----
        legends show residue number in 1-based index

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
    act = adict['multidihedral']
    act(command, traj, _top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype=dtype)

def calc_atomicfluct(traj=None, command="", *args, **kwd):
    dslist = DataSetList()
    dslist.set_py_free_mem(False)
    act = adict['atomicfluct']
    act(command, traj, dslist=dslist, *args, **kwd)
    # tag: print_output()
    act.print_output() # need to have this. check cpptraj's code
    return dslist[-1]

def calc_vector(traj=None, mask="", top=None, dtype='dataset', *args, **kwd): 
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
    from pytraj.actions.CpptrajActions import Action_Vector
    from pytraj.core.ActionList import ActionList

    dslist = DataSetList()
    _top = _get_top(traj, top)
    list_of_commands = _get_list_of_commands(mask)
    actlist = ActionList()

    for command in list_of_commands:
        act = Action_Vector()
        actlist.add_action(act, command, _top, dslist=dslist, *args, **kwd)
    actlist.do_actions(traj)

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
    if use_mass:
        command += " mass"

    act.read_input(command=command, top=_top, dslist=dslist)
    act.process(_top)

    for frame in _frame_iter_master(traj):
        # set Frame masses
        if use_mass:
            frame.set_frame_mass(_top)
        act.do_action(frame)
    return _get_data_from_dtype(dslist, dtype=dtype)

calc_COM = calc_center_of_mass = partial(_calc_vector_center, use_mass=True)

def calc_center_of_geometry(traj=None, command="", top=None, dtype='dataset'):
    _top = _get_top(traj, top)
    atom_mask_obj = _top(command)
    dslist = DataSetList()
    dslist.add_set("vector")
    #dslist.set_py_free_mem(False) # need this to avoid segmentation fault
    for frame in _frame_iter_master(traj):
        dslist[0].append(frame.center_of_geometry(atom_mask_obj))
    return _get_data_from_dtype(dslist, dtype=dtype)

calc_COG = calc_center_of_geometry

def calc_pairwise_rmsd(traj=None, command="", top=None, *args, **kwd):
    """return  DataSetList object
    Parameters
    ----------
    traj : Trajectory-like, iterable object
    command : mask (default=all atom) + extra command
        See `Notes` below for further info
    top : Topology, optional, default=None
    *args, **kwd: optional (dtype='dataset' | 'pyarray' | 'ndarray' | 'list')

    Examples:
        1. 
        # memory saving
        # * Create TrajectoryIterator to load only one Frame at a time
        traj = mdio.iterload("data/nogit/tip3p/md.trj", 
                            "./data/nogit/tip3p/tc5bwat.top")
        # we will load stripped-atom frames into memory only
        new_top = traj.top.strip_atoms("!@CA", copy=True)

        # passing `frame_iter` to `calc_pairwise_rmsd`
        # traj(0, 1000, mask='@CA') is equal to
        #     traj.frame_iter(start=0, stop=1000, mask='@CA')

        import pytraj.common_actions as pyca
        pyca.calc_pairwise_rmsd(traj(0, 1000, mask='@CA'), 
                                       top=new_top, dtype='ndarray')

        2. 
        # calculate pairwise rmsd for all frames using CA atoms
        dslist = calc_pairwise_rmsd(traj, "@CA")
        dslist.to_ndarray()
        dslist.tolist()

        3.
        # calculate pairwise rmsd for all frames using CA atoms, use `dme` (distance RMSD)
        # convert to numpy array
        arr_np = calc_pairwise_rmsd(traj, "@CA dme", dtype='ndarray')

        4.
        # calculate pairwise rmsd for all frames using CA atoms, nofit for RMSD
        # convert to numpy array
        arr_np = calc_pairwise_rmsd(traj, "@CA nofit", dtype='ndarray')

        5.
        # calculate pairwise rmsd for all frames using CA atoms
        # use symmetry-corrected RMSD, convert to numpy array
        arr_np = calc_pairwise_rmsd(traj, "@CA srmsd", dtype='ndarray')

    Notes
    -----
    * Same as `Analysis_Rms2d` in cpptraj (Amber15 manual, page 613)
    * Support `openmp`: require install `libcpptraj` with `openmp` flag
    * Memory: this calculation will make a copy of `traj`, it's better to use
        TrajectoryIterator
    * Command from cpptraj
        [<name>] [<mask>] [out <filename>]
        [dme | nofit | srmsd] [mass]
    * See full details
        from pytraj import info
        info("rms2d")
    """
    from pytraj.analyses.CpptrajAnalyses import Analysis_Rms2d
    if 'dtype' in kwd.keys():
        dtype = kwd['dtype']
        del kwd['dtype']
    else:
        dtype = None

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

    return _get_data_from_dtype(dslist, dtype)

def calc_density(traj=None, command="", top=None, *args, **kwd):
    # NOTE: trick cpptraj to write to file first and the reload
    if 'dtype' in kwd.keys():
        dtype = kwd['dtype']
        del kwd['dtype']
    else:
        dtype = None
    
    with goto_temp_folder():
        def _calc_density(traj, command, *args, **kwd):
            # TODO: update this method if cpptraj save data to DataSetList
            from pytraj.actions.CpptrajActions import Action_Density
        
            _top = _get_top(traj, top)
            dflist = DataFileList()
        
            tmp_filename = "tmp_pytraj_out.txt"
            command = "out " + tmp_filename + " " + command
            act = Action_Density()
            #with goto_temp_folder():
            act(command, traj, top=_top, dflist=dflist)
            act.print_output()
            dflist.write_all_datafiles()
            absolute_path_tmp = os.path.abspath(tmp_filename)
            return absolute_path_tmp

        dslist = DataSetList()
        fname = _calc_density(traj, command, *args, **kwd)
        dslist.read_data(fname)
        return _get_data_from_dtype(dslist, dtype)

def calc_temperatures(traj=None, command="", top=None):
    """return 1D python array of temperatures (from velocity) in traj
    if `frame` keyword is specified cpptraj/pytraj will take existing T

    Default = array of 0.0
    """
    from array import array as pyarray
    _top = _get_top(traj, top)
    dslist = calculate('temperature', traj, command, _top)
    return pyarray('d', dslist[0].tolist())

def calc_rmsd(traj=None, command="", ref=None, mass=False, 
              fit=True, top=None, dtype='pyarray',
              mode='pytraj'):
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
    dtype : data type, default='pyarray'
    mode : str {'pytraj', 'cpptraj'}
        if 'pytraj' a bit slower but original coords are not updated
            mask (command) can be string mask for atom index array
        if 'cpptraj': faster and coords for frame/traj are updated (rmsfit)
            only string mask for mask (command)
    Examples
    --------
    >>> from pytraj import io
    >>> from pytraj.common_actions import calc_rmsd
    >>> traj = io.load_sample_data("tz2")
    >>> calc_rmsd(traj, ":3-13@CA", ref=traj[0], mass=True, fit=True)
    >>> calc_rmsd(traj, ":3-13@CA", ref=traj[0], mass=True, fit=False)
    >>> calc_rmsd(traj, ":3-13@CA", ref=traj[0], mass=True, fit=False, mode='cpptraj')
    >>> calc_rmsd([traj, traj[-1]], ":3-13@CA", ref=traj[0], top=traj.top, mass=True, fit=False)

    """
    from array import array as pyarray
    from pytraj.datasets import DatasetDouble

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

    if mode == 'pytraj':
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
            ref.set_frame_mass(_top)
        _ref = Frame(ref, atm)
        for frame in _frame_iter_master(traj):
            if mass:
                # TODO : just need to set mass once
                frame.set_frame_mass(_top)
            if fit:
                _rmsd = frame.rmsd(ref, atommask=atm, use_mass=mass)
            else:
                _frame = Frame(frame, atm)
                _rmsd = _frame.rmsd_nofit(_ref, use_mass=mass)
            arr.append(_rmsd)
        if dtype == 'pyarray':
            return arr
        else:
            dset = DatasetDouble()
            dset.resize(len(arr))
            dset.values[:] = arr
            dset.legend = 'rmsd'
            return _get_data_from_dtype(dset, dtype=dtype)

    elif mode == 'cpptraj':
        if not isinstance(command, string_types):
            raise ValueError("only support string mask/command in mode=cpptraj")
        from pytraj.actions.CpptrajActions import Action_Rmsd
        act = Action_Rmsd()
        dslist = DataSetList()
        act(command, [ref, traj], top=_top, dslist=dslist)

        if dtype == 'pyarray':
            return pyarray('d', dslist[0].data)[1:]
        else:
            dset = DatasetDouble()
            dset.resize(dslist[0].size - 1)
            dset.values[:] = pyarray('d', dslist[0].data[1:])
            dset.legend = 'rmsd'
            return _get_data_from_dtype(dset, dtype=dtype)
    else:
        raise ValueError("mode = `pytraj` or `cpptraj`")


# alias for `calc_rmsd`
rmsd = calc_rmsd

def align_principal_axis(traj=None, command="*", top=None):
    # TODO : does not match with cpptraj output
    # rmsd_nofit ~ 0.5 for md1_prod.Tc5b.x, 1st frame
    """
    Notes
    -----
    apply for mutatble traj (Trajectory, Frame)
    """
    act = adict['principal']
    command += " dorotation"
    act(command, traj, top)

def closest(traj=None, command=None, top=None, *args, **kwd):
    """
    Parameters
    ----------
    traj : Trajectory-like | list of Trajectory-like/frames | frame_iter | chunk_iter
        traj could be anything as long as _frame_iter_master(traj) returns Frame
    command : str
        cpptraj command (below)
    top : Topology-like object, default=None, optional
    *args, **kwd: more arguments
        if dtype == 'dataset': return a tuple (new_traj, datasetlist)
            cpptraj only save data to DataSetList if `closestout` is specified (check example)
        if dtype != 'dataset': return only new_traj

    cpptraj command
    ---------------
    (from Amber15 manual, page 547: http://ambermd.org/doc12/Amber15.pdf)
    <# to keep> <mask> [noimage] [first | oxygen] [center]
    [closestout <filename>] [name <setname>] [outprefix <parmprefix>]
    
    where:
        <mask> Mask of atoms to search for closest waters around.
        [noimage] Do not perform imaging; only recommended if trajectory has previously
        been imaged.
        [first | oxygen] Calculate distances between all atoms in <mask> and the first atom
        of solvent only (recommended for standard water models as it will increase
        speed of calculation).
        [center] Search for waters closest to center of <mask> instead of each atom in
        <mask>.
        [closestout <filename>] Write information on the closest solvent molecules to
        <filename>.
        [outprefix <prefix>] Write corresponding topology to file with name prefix
        <prefix>.
        DataSet Aspects:
        [Frame] Frame number.
        [Mol] Original solvent molecule number.
        [Dist] Solvent molecule distance in Angstrom.
        [FirstAtm] First atom number of original solvent molecule.

    Similar to the strip command, but modify coordinate frame and topology by keeping only the specified number of
    closest solvent molecules to the region specified by the given mask. Solvent molecules can be determined
    automatically by cpptraj (by default residues named WAT, HOH, or TIP3).

    Examples
    --------
    >>> from pytraj import io
    >>> traj = io.load_sample_data ('tz2')
    >>> import pytraj.common_actions as pyca
    >>> # obtain new traj, keeping only closest 100 waters 
    >>> # to residues 1 to 13 (index starts from 1) by distance to the first atom of water
    >>> t = pyca.closest (traj, "100 :1-13 first")
    >>> # get new traj and get new DataSetList object to store more information
    >>> # (such as Frame number, original solvent molecule number, ...) (from cpptraj manual)
    >>> new_traj, dslist = pyca.closest (traj, "100 :1-13 first closestout test.out", dtype='dataset')
    >>> new_traj, dslist = pyca.closest (traj, "100 :1-13 first closestout test.out", dtype='dataframe')
    """

    from .actions.CpptrajActions import Action_Closest
    from pytraj.Trajectory import Trajectory
    from pytraj import DataSetList
    dslist = DataSetList()

    if 'dtype' in kwd.keys():
        dtype = kwd['dtype']
    else:
        dtype = None

    act = Action_Closest()
    fa = Trajectory()

    _top = _get_top(traj, top)
    new_top = Topology()
    new_top.py_free_mem = False # cpptraj will do
    if dtype and 'closestout' not in command:
        # trick cpptraj to dump data to DataSetList too
        command = command + " closestout tmp_pytraj_closestout.out"
    act.read_input(command, _top, dslist=dslist)
    act.process(_top, new_top)

    fa.top = new_top.copy()
    for frame in _frame_iter_master(traj):
        new_frame = Frame()
        new_frame.py_free_mem = False # cpptraj will do
        act.do_action(frame, new_frame)
        fa.append(new_frame.copy())

    if dtype:
        new_dslist = _get_data_from_dtype(dslist, dtype=dtype)
        return (fa, new_dslist)
    else:
        return fa

def native_contacts(traj=None, command="", top=None, dtype='dataset',
                    ref=None,
                    *args, **kwd):
    """
    Notes
    ----
    if `ref` is not None: first number in result corresponds to reference
    """
    from .actions.CpptrajActions import Action_NativeContacts
    act = Action_NativeContacts()
    dslist = DataSetList()

    _top = _get_top(traj, top)
    if ref is not None:
        act(command, [ref, traj], top=_top, dslist=dslist, *args, **kwd)
    else:
        act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype=dtype)

def nastruct(traj=None, command="", top=None, dtype='dataset',
                    *args, **kwd):
    """
    Examples
    --------
        dslist = nastruct(traj)
        dslist.groupby("major", mode='aspect') # information for major groove
        print (dslist.get_aspect())

    See Also
    --------
        Amber15 manual (http://ambermd.org/doc12/Amber15.pdf page 580)
    """
    # TODO: doc, rename method, move to seperate module?
    from .actions.CpptrajActions import Action_NAstruct
    act = Action_NAstruct()
    dslist = DataSetList()

    _top = _get_top(traj, top)
    act(command, traj, dslist=dslist, top=_top, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype=dtype)

def calc_grid(traj=None, command="", top=None, dtype='dataset',
                    *args, **kwd):
    """
    Examples
    --------
        dslist = calc_grid(traj)
    See Also
    --------
        Amber15 manual (http://ambermd.org/doc12/Amber15.pdf)
    """
    # TODO: doc, rename method, move to seperate module?
    from .actions.CpptrajActions import Action_Grid
    act = Action_Grid()
    dslist = DataSetList()

    # cpptraj require output
    command = "tmp_pytraj_grid_output.txt " + command
    _top = _get_top(traj, top)
    act(command, traj, dslist=dslist, top=_top, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype=dtype)

calc_grid.__doc__ = _long_manual.__grid__

def check_structure(traj=None, command="", top=None,
                    *args, **kwd):
    """
    Examples
    --------
    >>> check_structure(traj[0], top=traj.top)
    """
    from .actions.CpptrajActions import Action_CheckStructure
    act = Action_CheckStructure()

    # cpptraj require output
    _top = _get_top(traj, top)
    act(command, traj, top=_top, *args, **kwd)
