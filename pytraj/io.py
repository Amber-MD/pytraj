from __future__ import absolute_import
from .externals.six import string_types, PY3
from .Topology import Topology
from .TrajectoryIterator import TrajectoryIterator
from .data_sample.load_sample_data import load_sample_data
from .Frame import Frame
from .Trajectory import Trajectory
from .trajs.Trajin_Single import Trajin_Single
from .trajs.Trajout import Trajout
from .utils.check_and_assert import make_sure_exist, is_frame_iter
from .utils import goto_temp_folder
from .externals._load_HDF5 import load_hdf5
from .externals._pickle import to_pickle, read_pickle
from .externals._json import to_json, read_json
from .load_cpptraj_file import load_cpptraj_file
from ._shared_methods import _frame_iter_master
from .dataframe import to_dataframe
from ._set_silent import set_error_silent
from ._guess_filetype import _guess_filetype

try:
    from .externals._load_ParmEd import load_ParmEd, _load_chem
except:
    load_ParmEd = None

try:
    from pytraj.externals._load_pseudo_parm import load_pseudo_parm
except:
    load_pseudo_parm = None

# load mdtraj and MDAnalysis
from .externals._load_mdtraj import load_mdtraj 
from .externals._load_MDAnalysis import load_MDAnalysis

try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

__all__ = ['load', 'load_hdf5', 'write_traj', 'read_parm', 'write_parm', 'save']

EXTRA_LOAD_METHODS = {'HDF5' : load_hdf5, }

def load(*args, **kwd):
    """try loading and returning appropriate values"""

    if args and is_frame_iter(args[0]):
        return _load_from_frame_iter(*args, **kwd)

    if kwd:
        try:
            if is_frame_iter(kwd['filename']):
                return _load_from_frame_iter(*args, **kwd)
        except:
            pass

    if 'filename' in kwd.keys():
        make_sure_exist(kwd['filename'])
    else:
        make_sure_exist(args[0])

    if len(args) + len(kwd) == 1:
        if len(args) == 1:
            filename = args[0]
        else:
            filename = kwd[kwd.keys()[0]]
        # try to use cpptraj to load Topology
        top = readparm(*args, **kwd)
        if hasattr(top, 'is_empty') and top.is_empty():
            # load file
            # try loading
            try:
                filetype = _guess_filetype(filename) 
                new_object = EXTRA_LOAD_METHODS[filetype](*args, **kwd)
                return new_object
            except:
                raise ValueError("don't know how to load file/files")
        else:
            return top
    else:
        # load to Trajectory object
        return loadtraj(*args, **kwd)[:]

def _load_from_filelist(*args, **kwd):
    """return a list of Trajectory"""
    args_less = args[1:]
    if isinstance(args[0], (list, tuple)):
        mylist = args[0]
    elif isinstance(args[0], string_types):
        # "remd.x.*"
        from glob import glob
        mylist = sorted(glob(args[0]))
    else:
        raise ValueError()
    return [loadtraj(filename, *args_less, **kwd)[:] for filename in mylist]

def iterload(*args, **kwd):
    """return TrajectoryIterator object
    """
    if kwd and 'indices' in kwd.keys():
        raise ValueError("do not support indices for TrajectoryIterator loading")
    return loadtraj(*args, **kwd)

def _iterload_from_filelist(*args, **kwd):
    """return a list of TrajectoryIterator"""
    """return TrajectoryIterator object
    """
    args_less = args[1:]

    if kwd and 'indices' in kwd.keys():
        raise ValueError("do not support indices for TrajectoryIterator loading")

    if isinstance(args[0], (list, tuple)):
        mylist = args[0]
    elif isinstance(args[0], string_types):
        # "remd.x.*"
        from glob import glob
        mylist = sorted(glob(args[0]))
    else:
        raise ValueError()
    return [loadtraj(filename, *args_less, **kwd) for filename in mylist]

def loadtraj(filename=None, top=Topology(), indices=None):
    """load trajectory from filename
    Parameters
    ----------
    filename : str
    top : {str, Topology}
    indices : {None, list, array ...}

    Returns
    -------
    TrajectoryIterator : if indices is None
    or 
    Trajectory : if there is indices
    """
    if not isinstance(top, Topology):
        top = Topology(top)
    if top.is_empty():
        raise ValueError("can not load file without Topology or empty Topology")
    ts = TrajectoryIterator()
    ts.load(filename, top)

    if indices is not None:
        farray = Trajectory()
        farray.top = top.copy()
        for i in indices:
            farray.append(ts[i])
        return farray
    elif is_frame_iter(filename):
        return _load_from_frame_iter(filename, top)
    else:
        return ts

def _load_from_frame_iter(traj_frame_iter, top=None):
    if top is None or top.is_empty():
        raise ValueError("must provide non-empty Topology")
    fa = Trajectory(traj_frame_iter, top=top)
    return fa

def iterload_remd(filename, top=Topology(), T="300.0"):
    """Load remd trajectory for single temperature.
    Example: Suppose you have replica trajectoris remd.x.00{1-4}. 
    You want to load and extract only frames at 300 K, use this "load_remd" method

    Parameters
    ----------
    filename : str
    top : {str, Topology objecd}
    T : {int, float, str}, defaul="300.0"

    Returns
    ------
    TrajectoryIterator object
    """
    from pytraj import CpptrajState

    state = CpptrajState()
    # add keyword 'remdtraj' to trick cpptraj
    trajin = filename + ' remdtraj remdtrajtemp ' + str(T)
    state.toplist.add_parm(top)

    # load trajin, add "is_ensemble = False" to trick cpptraj
    # is_ensemble has 3 values: None, False and True
    state.add_trajin(trajin, is_ensemble=False)
    tlist = state.get_trajinlist()
    # get TrajectoryREMDIterator
    traj = tlist._getitem_remd(0)
    traj.top = state.toplist[0].copy()

    # use _tmpobj to hold CpptrajState(). If not, cpptraj will free memory
    traj._tmpobj = state
    return traj

def load_remd(filename, top=Topology(), T="300.0"):
    return iterload_remd(filename, top, T)[:]

def writetraj(filename="", traj=None, top=None, 
              fmt='UNKNOWN_TRAJ', indices=None,
              overwrite=False):
    """writetraj(filename="", traj=None, top=None, 
              ftm='UNKNOWN_TRAJ', indices=None):
    """
    # TODO : support list (tuple) of Trajectory, TrajectoryIterator or 
    # list of filenames
    #filename = filename.encode("UTF-8")

    if fmt == 'unknown':
        fmt = fmt.upper() + "_TRAJ"
    else:
        fmt = fmt.upper()

    if isinstance(top, string_types):
        top = Topology(top)

    if traj is None or top is None:
        raise ValueError("Need non-empty traj and top files")

    with Trajout(filename=filename, top=top, fmt=fmt, overwrite=overwrite) as trajout:
        if isinstance(traj, Frame):
            if indices is not None:
                raise ValueError("indices does not work with single Frame")
            trajout.writeframe(0, traj, top)
        else:
            if isinstance(traj, string_types):
                traj2 = load(traj, top)
            else:
                traj2 = traj

            if indices is None:
                # write all traj
                if isinstance(traj2, (Trajectory, TrajectoryIterator)):
                    for idx, frame in enumerate(traj2):
                        trajout.writeframe(idx, frame, top)
                elif isinstance(traj2, (list, tuple)):
                    # list, tuple
                    for traj3 in traj2:
                        for idx, frame in enumerate(traj3):
                            trajout.writeframe(idx, frame, top)
            else:
                if isinstance(traj2, (list, tuple)):
                    raise NotImplementedError("must be Trajectory or TrajectoryIterator instance")
                for idx in indices:
                    trajout.writeframe(idx, traj2[idx], top)


def writeparm(filename=None, top=None, fmt='AMBERPARM'):
    # TODO : add *args
    from pytraj.parms.ParmFile import ParmFile
    #filename = filename.encode("UTF-8")
    parm = ParmFile()
    parm.writeparm(filename=filename, top=top, fmt=fmt)

def readparm(filename):
    """return topology instance from reading filename"""
    #filename = filename.encode("UTF-8")
    set_error_silent(True)
    top = Topology(filename)
    set_error_silent(False)
    return top

def loadpdb_rcsb(pdbid):
    """load pdb file from rcsb website

    Parameters
    ----------
    pdbid : str

    Examples
    --------
        io.loadpdb_rcsb("2KOC") # popular RNA hairpin
    """

    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
    txt = urlopen(url).read()
    fname = "/tmp/tmppdb.pdb"
    with open(fname, 'w') as fh:
        if PY3:
            txt = txt.decode()
        fh.write(txt)
    traj = load(fname, fname)
    return traj

def load_single_frame(frame=None, top=None):
    """load single Frame"""
    return load(frame, top)[0]

def load_full_ParmEd(parmed_obj):
    """save and reload ParmEd object to pytraj object"""
    import os
    import tempfile

    name = "mytmptop"
    cwd = os.getcwd()
    directory_name = tempfile.mkdtemp()
    os.chdir(directory_name)
    parmed_obj.write_parm(name)
    top = load(name)
    os.remove(name)
    os.removedirs(directory_name)
    os.chdir(cwd)
    return top

# creat alias
write_traj = writetraj
save = writetraj
save_traj = writetraj
load_traj = loadtraj
read_parm = readparm
write_parm = writeparm
