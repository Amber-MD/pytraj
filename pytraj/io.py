from __future__ import absolute_import
from .externals.six import string_types, PY3
from .datafiles.load_sample_data import load_sample_data
from .utils.check_and_assert import ensure_exist, is_frame_iter
from .utils import goto_temp_folder
from .externals._load_HDF5 import load_hdf5
from .externals._pickle import to_pickle, read_pickle
from .externals._json import to_json, read_json
from .datasets.utils import load_datafile
from .datafiles.load_cpptraj_file import load_cpptraj_file
from ._shared_methods import _frame_iter_master
from ._set_silent import set_error_silent
from ._guess_filetype import _guess_filetype
from ._get_common_objects import _get_top
from .compat import zip

load_cpptraj_datafile = load_datafile

try:
    from .externals._load_ParmEd import load_ParmEd, _load_parmed
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

__all__ = ['load', 'iterload', 'load_remd', 'iterload_remd',
           '_load_from_filelist', '_iterload_from_filelist',
           'load_pdb_rcsb', 'load_pdb',
           'load_pseudo_parm', 'load_cpptraj_file',
           'load_datafile', 'load_hdf5',
           'load_sample_data',
           'load_ParmEd', 'load_full_ParmEd',
           'load_mdtraj',
           'load_MDAnalysis', 'load_MDAnalysisIterator',
           'load_topology', 'read_parm', 'write_parm', 
           'save', 'write_traj',
           'read_pickle', 'read_json',
           'to_pickle', 'to_json',
           ]

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
        ensure_exist(kwd['filename'])
    else:
        ensure_exist(args[0])

    if len(args) + len(kwd) == 1:
        if len(args) == 1:
            filename = args[0]
        else:
            filename = kwd[kwd.keys()[0]]
        # try to use cpptraj to load Topology
        top = read_parm(*args, **kwd)
        if hasattr(top, 'is_empty') and top.is_empty():
            try:
                # use ParmEd to load if cpptraj fails
                import parmed
                return load_pseudo_parm(parmed.load_file(args[0]))
            except:
                try:
                    # try to predict filetype and use proper loading method
                    filetype = _guess_filetype(filename) 
                    new_object = EXTRA_LOAD_METHODS[filetype](*args, **kwd)
                    return new_object
                except:
                    raise ValueError("don't know how to load file/files")
        else:
            return top
    else:
        # load to Trajectory object
        return load_traj(*args, **kwd)[:]

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
    return [load_traj(filename, *args_less, **kwd)[:] for filename in mylist]

def iterload(*args, **kwd):
    """return TrajectoryIterator object
    """
    if kwd and 'indices' in kwd.keys():
        raise ValueError("do not support indices for TrajectoryIterator loading")
    if kwd and 'engine' in kwd.keys() and kwd['engine'] == 'mdtraj':
        raise ValueError("do not support iterload with engine=='mdtraj'")
    return load_traj(*args, **kwd)

def _iterload_from_filelist(filename=None, top=None, force_load=False, *args, **kwd):
    """return a list of TrajectoryIterator"""

    if kwd and 'indices' in kwd.keys():
        raise ValueError("do not support indices for TrajectoryIterator loading")

    if isinstance(filename, (list, tuple)):
        trajnamelist = filename
    elif isinstance(filename, string_types):
        # "remd.x.*"
        from glob import glob
        trajnamelist = sorted(glob(filename))
    else:
        raise ValueError()

    if isinstance(top, (list, tuple)):
        toplist = top
    elif isinstance(top, string_types):
        # "remd.x.*"
        from glob import glob
        toplist = sorted(glob(top))
    else:
        raise ValueError()

    if len(trajnamelist) != len(toplist):
        if not force_load:
            raise ValueError("len of filename list is not equal to len of toplist")
        else:
            assert len(trajnamelist) > len(toplist), "toplist must have smaller len"
            last_top = toplist[-1]
            toplist += [last_top for _ in range(len(toplist), len(trajnamelist))] 

    return [load_traj(_filename, _top, *args, **kwd) 
            for _filename, _top  in zip(trajnamelist, toplist)]

def load_traj(filename=None, top=None, indices=None, engine='pytraj', *args, **kwd):
    """load trajectory from filename
    Parameters
    ----------
    filename : str
    top : {str, Topology}
    indices : {None, list, array ...}
    engine : str, {'pytraj', 'mdanalysis'}, default 'pytraj'
        if 'pytraj', use pytraj for iterload (return `TrajectoryIterator`)
        if 'mdanalysis', use this package (return `TrajectoryMDAnalysisIterator`)
    *args, **kwd: additional arguments, depending on `engine`

    Returns
    -------
    TrajectoryIterator : if indices is None and engine='pytraj'
    or 
    Trajectory : if there is indices
    or TrajectoryMDAnalysisIterator if engine='mdanalysis'
    """
    if 'frame_slice' in kwd.keys() and not engine == 'pytraj':
        raise KeyError("only support frame_slice in engine mode = 'pytraj'")

    engine = engine.lower()

    if engine == 'pytraj':
        from .Topology import Topology
        from .TrajectoryIterator import TrajectoryIterator
        from .Trajectory import Trajectory

        if not isinstance(top, Topology):
            top = Topology(top)
        if top.is_empty():
            raise ValueError("can not load file without Topology or empty Topology")
        ts = TrajectoryIterator(top=top)

        if 'frame_slice' in kwd.keys():
            ts.load(filename, frame_slice=kwd['frame_slice'])
        else:
            ts.load(filename)

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
    elif engine == 'mdanalysis':
        from MDAnalysis import Universe as U
        if top is None:
            top = filename
        return load_MDAnalysisIterator(U(top, filename, *args, **kwd))
    elif engine == 'mdtraj':
        import mdtraj as md
        if top is None:
            top = filename
        return load_mdtraj(md.load(filename, top=top, *args, **kwd))
    else:
        raise NotImplementedError("support only {'pytraj', 'mdanlaysis', 'mdtraj'} engines")


def _load_from_frame_iter(traj_frame_iter, top=None):
    from .Trajectory import Trajectory
    if top is None or top.is_empty():
        if hasattr(traj_frame_iter, 'top'):
            top = traj_frame_iter.top
        else:
            raise ValueError("must provide non-empty Topology")
    fa = Trajectory(traj_frame_iter, top=top)
    return fa

def iterload_remd(filename, top=None, T="300.0"):
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

def load_remd(filename, top=None, T="300.0"):
    return iterload_remd(filename, top, T)[:]

def write_traj(filename="", traj=None, top=None, 
              format='unknown_traj', indices=None,
              overwrite=False, more_args="", 
              *args, **kwd):
    """write Trajectory-like, list of trajs, frames, ... to file/files

    Suppot file extensions
    ----------------------
    .crd, .nc, .rst7, .ncrst, .dcd, .pdb, .mol2, .binpos, .trr, .sqm
    if extension or format is not specify correctly, 
    cpptraj will use Amber Trajectory format (.crd)

    Examples
    --------
    >>> from pytraj import io
    >>> traj = io.load_sample_data()
    >>> io.write_traj("t.nc", traj) # write to amber netcdf file
    >>> # write to multi pdb files (t.pdb.1, t.pdb.2, ...)
    >>> io.write_traj("t.pdb", traj, overwrite=True, more_args='multi')
    >>> # write all frames to single pdb file and each frame is seperated by "MODEL" word
    >>> io.write_traj("t.pdb", traj, overwrite=True, more_args='model')
    >>> # write to DCD file
    >>> io.write_traj("test.dcd", traj)
    >>> # set nobox for trajout
    >>> io.write_traj("test.nc", traj, more_args='nobox')

    See Also
    --------
    Amber15 manual (http://ambermd.org/doc12/Amber15.pdf, page 542)

    Excerpt from this manual
    -----------------------
    Options for pdb format:
    [model | multi] [dumpq | parse | vdw] [chainid <ID>]
    [pdbres] [pdbatom] [pdbv3] [teradvance]

    Options for mol2 format:
    [single | multi]

    Options for SQM input format:
    [charge <c>]
    """
    from .Frame import Frame
    from .trajs.Trajout import Trajout

    if format.upper() == 'UNKNOWN':
        format= format.upper() + "_TRAJ"
    else:
        format= format.upper()

    _top = _get_top(traj, top)
    if _top is None:
        raise ValueError("must provide Topology")

    if traj is None or _top is None:
        raise ValueError("Need non-empty traj and top files")

    with Trajout(filename=filename, top=_top, format=format, 
                 overwrite=overwrite, more_args=more_args,
                 *args, **kwd) as trajout:
        if isinstance(traj, Frame):
            if indices is not None:
                raise ValueError("indices does not work with single Frame")
            trajout.writeframe(0, traj, _top)
        else:
            if isinstance(traj, string_types):
                traj2 = load(traj, _top)
            else:
                traj2 = traj

            if indices is not None:
                if isinstance(traj2, (list, tuple)):
                    raise NotImplementedError("must be Trajectory or TrajectoryIterator instance")
                for idx in indices:
                    trajout.writeframe(idx, traj2[idx], _top)

            else:
                for idx, frame in enumerate(_frame_iter_master(traj2)):
                    trajout.writeframe(idx, frame, _top)


def write_parm(filename=None, top=None, format='AMBERPARM'):
    # TODO : add *args
    from pytraj.parms.ParmFile import ParmFile
    #filename = filename.encode("UTF-8")
    parm = ParmFile()
    parm.writeparm(filename=filename, top=top, format=format)

def read_parm(filename):
    from .Topology import Topology
    """return topology instance from reading filename"""
    #filename = filename.encode("UTF-8")
    set_error_silent(True)
    top = Topology(filename)
    set_error_silent(False)
    return top

# creat alias
load_topology = read_parm

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

def download_PDB(pdbid, location="./", overwrite=False):
    """download pdb to local disk

    Return
    ------
    None

    Notes
    -----
    this method is different from `parmed.download_PDB`, which return a `Structure` object
    """
    import os
    fname = location + pdbid + ".pdb"
    if os.path.exists(fname) and not overwrite:
        raise ValueError("must set overwrite to True")

    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
    txt = urlopen(url).read()
    with open(fname, 'w') as fh:
        if PY3:
            txt = txt.decode()
        fh.write(txt)

# create alias
load_pdb_rcsb = loadpdb_rcsb

def load_pdb(pdb_file):
    """return a Trajectory object"""
    return load_traj(pdb_file, pdb_file)

def load_single_frame(frame=None, top=None, index=0):
    """load single Frame"""
    return load(frame, top)[index]

def load_full_ParmEd(parmed_obj):
    """save and reload ParmEd object to pytraj object"""
    import os
    import tempfile

    with goto_temp_folder():
        name = "mytmptop"
        if hasattr(parmed_obj, 'write_parm'):
            parmed_obj.write_parm(name)
        elif hasattr(parmed_obj, 'write_pdb'):
            name = name + ".pdb"
            parmed_obj.write_pdb(name)
        top = load_topology(name)
    return top

def load_MDAnalysisIterator(u):
    from .trajs.TrajectoryMDAnalysisIterator import TrajectoryMDAnalysisIterator
    return TrajectoryMDAnalysisIterator(u)

# creat alias
save = write_traj
save_traj = write_traj

def get_coordinates(an_object, top=None):
    '''return 3D-ndarray coordinates of `an_object`
    Parameters
    ----------
    an_object : could be anything having Frame info
        a Trajectory, TrajectoryIterator,
        a frame_iter, FrameIter, ...
    top : optional Topology if `an_object` does not have this information

        This method is designed to load coordinates with minimum memory requirement
    '''
    if hasattr(an_object, 'xyz'):
        return an_object.xyz[:]
    elif is_frame_iter(an_object):
        return _load_from_frame_iter(an_object, top=top).xyz[:]
