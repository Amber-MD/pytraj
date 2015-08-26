from __future__ import absolute_import
from .externals.six import string_types, PY3
from .datafiles.load_sample_data import load_sample_data
from .utils.check_and_assert import ensure_exist, is_frame_iter
from .utils import goto_temp_folder
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

# load mdtraj and MDAnalysis
from .externals._load_mdtraj import load_mdtraj
from .externals._load_MDAnalysis import load_MDAnalysis

try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

__all__ = ['load',
           'iterload',
           'load_remd',
           'iterload_remd',
           '_load_from_filelist',
           '_iterload_from_filelist',
           'load_pdb_rcsb',
           'load_pdb',
           'load_cpptraj_file',
           'load_datafile',
           'load_hdf5',
           'load_sample_data',
           'load_ParmEd',
           'load_mdtraj',
           'load_MDAnalysis',
           'load_MDAnalysisIterator',
           'load_topology',
           'read_parm',
           'write_parm',
           'save',
           'write_traj',
           'read_pickle',
           'read_json',
           'to_pickle',
           'to_json', ]


def load(*args, **kwd):
    """try loading and returning appropriate values

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load('traj.nc', '2koc.parm7')
    >>> traj = pt.load('traj.mol2')
    >>> traj = pt.load('traj.pdb')

    Notes
    -----
    load(filename, top, use_numpy=True) will return `pytraj.api.Trajectory`
    """

    if args and is_frame_iter(args[0]):
        return _load_from_frame_iter(*args, **kwd)

    if kwd:
        try:
            if is_frame_iter(kwd['filename']):
                return _load_from_frame_iter(*args, **kwd)
        except:
            pass

    if 'filename' in kwd.keys():
        filename = kwd['filename']
    else:
        filename = args[0]

    if filename.startswith('http://') or filename.startswith('https://'):
        return load_ParmEd(filename, as_traj=True, structure=True)
    else:
        ensure_exist(filename)
        # load to TrajectoryIterator object first
        traj = load_traj(*args, **kwd)
        if 'use_numpy' in kwd.keys() and kwd['use_numpy']:
            from pytraj.api import Trajectory
            return Trajectory(xyz=traj.xyz, top=traj.top)
        else:
            return traj[:]


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

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.iterload('traj.nc', '2koc.parm7')
    >>> traj = pt.iterload(['traj0.nc', 'traj1.nc'], '2koc.parm7')
    >>> traj = pt.iterload('./traj*.nc', '2koc.parm7')
    """
    if kwd and 'indices' in kwd.keys():
        raise ValueError(
            "do not support indices for TrajectoryIterator loading")
    if kwd and 'engine' in kwd.keys() and kwd['engine'] == 'mdtraj':
        raise ValueError("do not support iterload with engine=='mdtraj'")
    return load_traj(*args, **kwd)


def _load_netcdf(filename, top, indices=None, engine='scipy'):
    from pytraj import api
    traj = api.Trajectory(top=top)

    if engine == 'scipy':
        from scipy import io
        fh = io.netcdf_file(filename, mmap=False)
        data = fh.variables['coordinates'].data
        clen = fh.variables['cell_lengths'].data
        cangle = fh.variables['cell_angles'].data
    if engine == 'netcdf4':
        import netCDF4
        fh = netCDF4.Dataset(filename)
        data = fh.variables['coordinates']
        clen = fh.variables['cell_lengths']
        cangle = fh.variables['cell_angles']
    if indices is None:
        traj.xyz = data
    else:
        traj.xyz = data[indices]
    if traj.xyz.itemsize != 8:
        traj.xyz = traj.xyz.astype('f8')
    traj._append_unitcells((clen, cangle))
    return traj


def _iterload_from_filelist(filename=None,
                            top=None,
                            force_load=False, *args, **kwd):
    """return a list of TrajectoryIterator"""

    if kwd and 'indices' in kwd.keys():
        raise ValueError(
            "do not support indices for TrajectoryIterator loading")

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
            raise ValueError(
                "len of filename list is not equal to len of toplist")
        else:
            assert len(trajnamelist) > len(
                toplist), "toplist must have smaller len"
            last_top = toplist[-1]
            toplist += [
                last_top for _ in range(len(toplist), len(trajnamelist))
            ]

    return [load_traj(_filename, _top, *args, **kwd)
            for _filename, _top in zip(trajnamelist, toplist)]


def load_traj(filename=None,
              top=None,
              indices=None,
              engine='pytraj', *args, **kwd):
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
        from .api import Trajectory

        if not isinstance(top, Topology):
            top = Topology(top)
        if top.is_empty():
            top = Topology(filename)
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
        raise NotImplementedError(
            "support only {'pytraj', 'mdanlaysis', 'mdtraj'} engines")


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
    from pytraj import Trajectory
    itertraj = iterload_remd(filename, top, T)
    return Trajectory(itertraj, top=itertraj.top)


def write_traj(filename="",
               traj=None,
               top=None,
               format=None,
               indices=None,
               overwrite=False,
               mode="", *args, **kwd):
    """write Trajectory-like, list of trajs, frames, ... to file/files

    Suppot file extensions
    ----------------------
    .crd, .nc, .rst7, .ncrst, .dcd, .pdb, .mol2, .binpos, .trr, .sqm
    if extension or format is not specify correctly, 
    cpptraj will use Amber Trajectory format (.crd)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = io.load_sample_data()
    >>> pt.write_traj("t.nc", traj, overwrite=True) # write to amber netcdf file

    >>> # write to multi pdb files (t.pdb.1, t.pdb.2, ...)
    >>> pt.write_traj("t.pdb", traj, overwrite=True, mode='multi')

    >>> # write all frames to single pdb file and each frame is seperated by "MODEL" word
    >>> pt.write_traj("t.pdb", traj, overwrite=True, mode='model')

    >>> # write to DCD file
    >>> pt.write_traj("test.dcd", traj, overwrite=True)
    """
    from .Frame import Frame
    from .trajs.Trajout import Trajout

    if format in [None, '']:
        # use cpptraj default format (amber)
        format = 'unknown_traj'
    if format.upper() == 'UNKNOWN':
        format = format.upper() + "_TRAJ"
    else:
        format = format.upper()

    _top = _get_top(traj, top)
    if _top is None:
        raise ValueError("must provide Topology")

    if traj is None or _top is None:
        raise ValueError("Need non-empty traj and top files")

    with Trajout(filename=filename,
                 top=_top,
                 format=format,
                 overwrite=overwrite,
                 mode=mode, *args, **kwd) as trajout:
        if isinstance(traj, Frame):
            if indices is not None:
                raise ValueError("indices does not work with single Frame")
            trajout.write(0, traj, _top)
        else:
            if isinstance(traj, string_types):
                traj2 = iterload(traj, _top)
            else:
                traj2 = traj

            if indices is not None:
                if isinstance(traj2, (list, tuple)):
                    raise NotImplementedError(
                        "must be Trajectory or TrajectoryIterator instance")
                for idx in indices:
                    trajout.write(idx, traj2[idx], _top)

            else:
                for idx, frame in enumerate(_frame_iter_master(traj2)):
                    trajout.write(idx, frame, _top)


def write_parm(filename=None, top=None, format='AMBERPARM'):
    # TODO : add *args
    from pytraj.parms.ParmFile import ParmFile
    #filename = filename.encode("UTF-8")
    parm = ParmFile()
    parm.writeparm(filename=filename, top=top, format=format)


def load_topology(filename, **kwd):
    """
    load Topology from filename or from url
    >>> import pytraj as pt
    >>> pt.load_topology("./data/tz2.ortho.parm7")
    >>> pt.load_topology("http://ambermd.org/tutorials/advanced/tutorial1/files/polyAT.pdb")
    """
    if filename.startswith('http://') or filename.startswith('https://'):
        return _load_url(filename)
    else:
        from .Topology import Topology
        """return topology instance from reading filename"""
        #filename = filename.encode("UTF-8")
        set_error_silent(True)
        top = Topology(filename)
        set_error_silent(False)
        return top

# creat alias
read_parm = load_topology


def _load_url(url):
    """load Topology from url
    """
    from .Topology import Topology

    txt = urlopen(url).read()
    fname = "/tmp/tmppdb.pdb"
    with open(fname, 'w') as fh:
        if PY3:
            txt = txt.decode()
        fh.write(txt)
    return Topology(fname)


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
    """load a single Frame"""
    return iterload(frame, top)[index]


load_frame = load_single_frame


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
