"""this file has commonly used actions such as rmsd calculation, 
randomizeions, strip atoms, ..."""

from __future__ import print_function, absolute_import
import os
from glob import glob
from pytraj.Topology import Topology
from .core.TopologyList import TopologyList
from .ArgList import ArgList
from pytraj.Frame import Frame
from pytraj.Trajectory import Trajectory
from pytraj.DataSetList import DataSetList
from pytraj._shared_methods import _frame_iter as frame_iter
from pytraj._set_silent import set_world_silent
from pytraj.compat import set

# external
from pytraj.externals.six import string_types

try:
    from pytraj.externals.magic import from_file as file_type_info
except ImportError:
    file_type_info = None

__all__ = ['to_amber_mask', 'from_legends_to_indices', 'info',
           'get_atts', 'merge_trajs']

def to_amber_mask(txt, mode=None):
    import re
    """Convert something like 'ASP_16@OD1-ARG_18@N-H to ':16@OD1 :18@H'

    Parameters
    ----------
    txt : str | list/tuple of string | array-like of integer
    mode : str, default=None
        if mode='int_to_str': convert integer array to Amber mask
            (good for converting indices to atom mask string to be used with cpptraj)

    Examples
    --------
        to_amber_mask('ASP_16@OD1-ARG_18@N-H') # get ':16@OD1 :18@H'
        to_amber_mask(range(0, 10, 3), mode='int_to_str') # return `@1,4,7`
    """

    if mode is None:
        if isinstance(txt, string_types):
            txt = txt.replace("_", ":")
            return " ".join(re.findall(r"(:\d+@\w+)", txt))
        elif isinstance(txt, (list, tuple)):
            # list is mutable
            txt_copied = txt[:]
            for i, _txt in enumerate(txt):
                txt_copied[i] = to_amber_mask(_txt)
            return txt_copied
        else:
            raise NotImplementedError()
    elif mode == 'int_to_str':
        # need to add +1 since cpptraj's mask uses starting index of 1
        my_long_str = ",".join(str(i+1) for i in txt)
        return "@" + my_long_str
    else:
        raise NotImplementedError()

def from_legends_to_indices(legends, top):
    """return somethine like "ASP_16@OD1-ARG_18@N-H" to list of indices

    Parameters
    ----------
    legends : str
    top : Topology
    """
    mask_list = to_amber_mask(legends)
    index_list = []
    for m in mask_list:
        index_list.append(top(m).indices)
    return index_list

def info(obj=None):
    """get `help` for obj
    Useful for Actions and Analyses
    
    Since we use `set_worl_silent` to turn-off cpptraj' stdout, we need 
    to turn on to use cpptraj's help methods
    """
    from pytraj import adict, analdict
    adict_keys = adict.keys()
    anal_keys = analdict.keys()

    if obj is None:
        print ("action's keys", adict_keys)
        print ("analysis' keys", anal_keys)
    else:
        if isinstance(obj, string_types):
            if obj in adict.keys():
                # make Action object
                _obj = adict[obj]
            elif obj in analdict.keys():
                # make Analysis object
                _obj = analdict[obj]
            else:
                raise ValueError("keyword must be an Action or Analysis")
        else:
            # assume `obj` hasattr `help`
            _obj = obj

        if hasattr(_obj, 'help'):
            set_world_silent(False)
            _obj.help()
            set_world_silent(True)
        elif hasattr(_obj, 'info'):
            set_world_silent(False)
            _obj.info()
            set_world_silent(True)
        elif 'calc_' in _obj.__name__:
            key = _obj.__name__.split("_")[-1]
            set_world_silent(False)
            adict[key].help()
            set_world_silent(True)
        elif hasattr(_obj, '__doc__'):
            print (_obj.__doc_)
        else:
            raise ValueError("object does not have `help` method")

def get_action_dict():
    from pytraj.actions import allactions
    actdict = {}
    for key in allactions.__dict__.keys():
        if "Action_" in key:
            act = key.split("Action_")[1]
            # add Action classes
            actdict[act] = allactions.__dict__["Action_" + act]
    return actdict

# add action_dict
action_dict = get_action_dict()

def show_code(func, get_txt=False):
    """show code of func or module"""
    import inspect
    txt = inspect.getsource(func)
    if not get_txt:
        print (txt)
    else:
        return txt

def get_atts(obj):
    """get methods and atts from obj but excluding special methods __"""
    atts_dict = dir(obj)
    return [a for a in atts_dict if not a.startswith("__")]

def merge_trajs(traj1, traj2, start_new_mol=True, n_frames=None):
    """
    
    Examples
    --------
       # from two Trajectory or TrajectoryIterator
       traj3 = merge_trajs(traj1, traj2)
       assert traj3.n_frames == traj1.n_frames == traj2.n_frames
       assert traj3.n_atoms == traj1.n_atoms + traj2.n_atoms
       import numpy as np
       assert np.any(traj3.xyz, np.vstack(tra1.xyz,  traj2.xyz)) == True

       # from frame_iter for saving memory
       traj3 = merge_trajs((traj1(0, 10, 2), traj1.top), 
                           (traj2(100, 110, 2), traj2.top), n_frames=6)

    Notes
    -----
    Code might be changed
    """
    from pytraj.compat import zip
    import numpy as np

    if isinstance(traj1, (list, tuple)):
        n_frames_1 = n_frames
        top1 = traj1[1]
        _traj1 = traj1[0]
    else:
        n_frames_1 = traj1.n_frames
        top1 = traj1.top
        _traj1 = traj1

    if isinstance(traj2, (list, tuple)):
        n_frames_2 = n_frames
        top2 = traj2[1] # example: (traj(0, 5), traj.top)
        _traj2 = traj2[0]
    else:
        n_frames_2 = traj2.n_frames
        top2 = traj2.top
        _traj2 = traj2

    if n_frames_1 != n_frames_2:
        raise ValueError("must have the same n_frames")

    if top1.n_atoms != top2.n_atoms:
        raise ValueError("must have the same n_atoms")

    traj = Trajectory()
    traj._allocate(n_frames_1, top1.n_atoms + top2.n_atoms)

    # merge Topology
    top = top1.copy()
    if start_new_mol:
        top.start_new_mol()
    top.join(top2)
    traj.top = top

    # update coords
    for f1, f2, frame in zip(_traj1, _traj2, traj):
        frame.xyz = np.vstack((f1.xyz, f2.xyz))

    return traj

def find_libcpptraj(**kwd):
    return find_library('cpptraj', **kwd)

def find_library(libname, unique=False):
    """return a list of all library files"""
    paths = os.environ.get('LD_LIBRARY_PATH', '').split(':')
    lib_path_list = []
    key = "lib" + libname + "*"

    for path in paths:
        path = path.strip()
        fnamelist = glob(os.path.join(path, key))
        for fname in fnamelist: 
            if os.path.isfile(fname):
                lib_path_list.append(fname)

    if not lib_path_list:
        return None
    else:
        if unique:
            return set(lib_path_list)
        else:
            return lib_path_list
