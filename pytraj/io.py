from __future__ import absolute_import
import os
import numpy as np

from .externals.six import string_types, PY3
from .datafiles.load_samples import load_sample_data
from .externals._pickle import to_pickle, read_pickle
from .externals._json import to_json, read_json
from .datafiles.load_cpptraj_file import load_cpptraj_file
from ._shared_methods import iterframe_master
from .cyutils import _fast_iterptr as iterframe_from_array
from .cpp_options import set_error_silent
from ._get_common_objects import _get_topology
from .topology import Topology, ParmFile
from .trajectory import Trajectory
from .trajectory_iterator import TrajectoryIterator

from .externals._load_ParmEd import load_ParmEd

from .decorators import ensure_exist

try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

__all__ = ['load',
           'iterload',
           'load_remd',
           'iterload_remd',
           'load_pdb_rcsb',
           'load_cpptraj_file',
           'load_sample_data',
           'load_ParmEd',
           'load_topology',
           'read_parm',
           'write_parm',
           'save',
           'write_traj',
           'read_pickle',
           'read_json',
           'to_pickle',
           'to_json', ]


def load(filename, top=None, frame_indices=None, mask=None):
    """try loading and returning appropriate values. See example below.

    Parameters
    ----------
    filename : str, Trajectory filename
    top : Topology filename or a Topology
    frame_indices : {None, array_like}, default None
        only load frames with given number given in frame_indices
    mask : {str, None}, default None
        if None: load coordinates for all atoms
        if string, load coordinates for given atom mask 

    Returns
    -------
    pytraj.Trajectory

    Notes
    -----
    - For further slicing options, see pytraj.TrajectoryIterator (created by ``pytraj.iterload``)

    - Also see `pytraj.iterload` for loading a series of trajectories that don't fit to
      memory

    See also
    --------
    iterload

    Examples
    --------
    >>> import pytraj as pt
    >>> # load netcdf file with given amber parm file
    >>> traj = pt.load('traj.nc', '2koc.parm7') # doctest: +SKIP

    >>> # load netcdf file with given amber parm file
    >>> traj = pt.load('traj.nc', '2koc.parm7') # doctest: +SKIP

    >>> # load mol2 file
    >>> traj = pt.load('traj.mol2') # # doctest: +SKIP

    >>> # load pdb file
    >>> traj = pt.load('traj.pdb') # doctest: +SKIP

    >>> # load given frame numbers
    >>> traj = pt.load('traj.nc', top='2koc.parm7', frame_indices=[0, 3, 5, 12, 20]) # doctest: +SKIP
    >>> traj = pt.load('traj.nc', top='2koc.parm7', frame_indices=[0, 3, 5, 12, 20], mask='!@H=') # doctest: +SKIP

    >>> # load with frame slice
    >>> traj = pt.load('traj.nc', top='2koc.parm7', frame_indices=slice(0, 10, 2)) # doctest: +SKIP
    >>> # which is equal to:
    >>> traj = pt.load('traj.nc', top='2koc.parm7', frame_indices=range(0, 10, 2)) # doctest: +SKIP
    """
    if isinstance(filename, string_types) and filename.startswith(
            'http://') or filename.startswith('https://'):
        return load_ParmEd(filename, as_traj=True, structure=True)
    else:
        # load to TrajectoryIterator object first
        traj = load_traj(filename, top)

        # do the slicing and other thinkgs if needed.
        if isinstance(frame_indices, tuple):
            frame_indices = list(frame_indices)
        if frame_indices is None and mask is None:
            # load all
            return traj[:]
        elif frame_indices is None and mask is not None:
            # load all frames with given mask
            # eg: traj['@CA']
            return traj[mask]
        elif frame_indices is not None and mask is None:
            # eg. traj[[0, 3, 7]]
            return traj[frame_indices]
        else:
            # eg. traj[[0, 3, 7], '@CA']
            return traj[frame_indices, mask]


def iterload(*args, **kwd):
    """return TrajectoryIterator object

    Parameters
    ----------
    filename: {str, list-like of filenames, pattern}
        input trajectory filename(s). You can use a single filename, a list of filenames
        or a pattern.
    top : {str, Topology}
        input topology. If str, pytraj will load from disk to Topology first
    frame_slice: tuple or list of tuple
        specify start, stop, step for each trajectory you want to read.

        cpptraj input::
            
            trajin traj0.nc 1 10
            trajin traj1.nc

        In pytraj, corresponding frame_slice=[(0, 10), (0, -1)]

    Returns
    -------
    pytraj.TrajectoryIterator

    Examples
    --------
    >>> import pytraj as pt
    >>> # load a single filename
    >>> traj = pt.iterload('traj.nc', '2koc.parm7') # doctest: +SKIP

    >>> # load from a list of filenames
    >>> traj = pt.iterload(['traj0.nc', 'traj1.nc'], '2koc.parm7') # doctest: +SKIP

    >>> # load all files with a given pattern (sorted)
    >>> traj = pt.iterload('./traj*.nc', '2koc.parm7') # doctest: +SKIP

    >>> # load from a list of files with given frame step
    >>> # for each file, only register to load from frame 0 to 9 (skip 10), skip every 2 frames
    >>> traj = pt.iterload(['traj0.nc', 'traj1.nc'], '2koc.parm7', frame_slice=[(0, 10, 2),]*2) # doctest: +SKIP

    >>> # load from frame 0 to 9 for `traj0.nc`
    >>> # load all frames from `traj1.nc`
    >>> traj = pt.iterload(['traj0.nc', 'traj1.nc'], '2koc.parm7', frame_slice=[(0, 10), (0, -1)] # doctest: +SKIP

    Notes
    -----
    Unlike `pytraj.load`, you can not arbitarily set `frame_indices`. If you want to do
    so, first load trajectories to TrajectoryIterator object, then do fancy slicing

    >>> import pytraj as pt
    >>> # register to load traj.nc from 0-th to 99-th frame
    >>> traj = pt.iterload('traj.nc', 'prmtop', frame_slice=(0, 100)]) # doctest: +SKIP
    >>> # do fancy indexing to load specific frames to memory
    >>> traj[[0, 8, 3, 50, 7]] # doctest: +SKIP

    >>> # load to disk with given mask
    >>> traj[[0, 8, 3, 50, 7], '!@H='] # doctest: +SKIP
    """
    if kwd and 'frame_indices' in kwd.keys():
        raise ValueError(
            "do not support indices for TrajectoryIterator loading")
    return load_traj(*args, **kwd)


def _load_netcdf(filename, top,
                 frame_indices=None,
                 engine='scipy'):  # pragma: no cover
    '''simply read all data to memory. Use this if you want to load data few times
    faster (and  you know what you are doing).
    '''
    from pytraj import trajectory
    traj = trajectory.Trajectory(top=top)

    if engine == 'scipy':
        from scipy import io
        fh = io.netcdf_file(filename, mmap=False)
        data = fh.variables['coordinates'].data
        try:
            clen = fh.variables['cell_lengths'].data
            cangle = fh.variables['cell_angles'].data
        except KeyError:
            clen = None
            cangle = None
    if engine == 'netcdf4':
        import netCDF4
        fh = netCDF4.Dataset(filename)
        data = fh.variables['coordinates']
        try:
            clen = fh.variables['cell_lengths']
            cangle = fh.variables['cell_angles']
        except KeyError:
            clen = None
            cangle = None
    if clen is not None and cangle is not None:
        clen = np.ascontiguousarray(clen, dtype='f8')
        cangle = np.ascontiguousarray(cangle, dtype='f8')
    if frame_indices is None:
        traj.xyz = np.ascontiguousarray(data, dtype='f8')
    else:
        traj.xyz = np.ascontiguousarray(data, dtype='f8')[frame_indices]
        if clen is not None:
            clen = clen[frame_indices]
            cangle = cangle[frame_indices]
    if clen is not None:
        traj._append_unitcells((clen, cangle))
    return traj


def load_traj(filename=None, top=None, *args, **kwd):
    """load trajectory from filename

    Parameters
    ----------
    filename : str
    top : {str, Topology}
    frame_indices : {None, list, array ...}
    *args, **kwd: additional arguments, depending on `engine`

    Returns
    -------
    TrajectoryIterator : if frame_indices is None
    Trajectory : if there is indices
    """
    if isinstance(top, string_types):
        top = load_topology(top)
    if top is None or top.is_empty():
        top = load_topology(filename)
    ts = TrajectoryIterator(top=top)

    if 'frame_slice' in kwd.keys():
        ts.load(filename, frame_slice=kwd['frame_slice'])
    else:
        ts.load(filename)

    return ts


def _load_from_frame_iter(iterable, top=None):
    '''
    '''
    return Trajectory.from_iterable(iterable, top)


def iterload_remd(filename, top=None, T="300.0"):
    """Load temperature remd trajectory for single temperature.
    Example: Suppose you have replica trajectoris remd.x.00{1-4}. 
    You want to load and extract only frames at 300 K, use this method

    Parameters
    ----------
    filename : str
    top : {str, Topology}
    T : {float, str}, default=300.0

    Returns
    -------
    pytraj.traj.TrajectoryCpptraj

    Notes
    -----

    """
    from pytraj.core.cpp_core import CpptrajState, Command
    dispatch = Command.dispatch

    state = CpptrajState()

    # add keyword 'remdtraj' to trick cpptraj
    trajin = ' '.join(('trajin', filename, 'remdtraj remdtrajtemp', str(T)))
    if isinstance(top, string_types):
        top = read_parm(top)
    else:
        top = top
    state.data.add_set('topology', 'remdtop')
    # set topology
    state.data['remdtop']._top = top

    # load trajin
    dispatch(state, trajin)
    dispatch(state, 'loadtraj name remdtraj')

    # state.data.remove_set(state.data['remdtop'])
    traj = state.data[-1]

    # assign state to traj to avoid memory free
    traj._base = state
    return traj


def load_remd(filename, top=None, T="300.0"):
    return iterload_remd(filename, top, T)[:]


def write_traj(filename="",
               traj=None,
               top=None,
               frame_indices=None,
               overwrite=False,
               options=""):
    """write Trajectory-like or iterable object to trajectory file

    Parameters
    ----------
    filename : str
    traj : Trajectory-like or iterator that produces Frame or 3D ndarray with shape=(n_frames, n_atoms, 3)
    top : Topology, optional, default: None
    frame_indices: array-like or iterator that produces integer, default: None
        If not None, only write output for given frame indices
    overwrite: bool, default: False
    options : str, additional keywords for extention='.pdb'. See examples.
        
    Notes
    -----
    ===================  =========
    Format               Extension
    ===================  =========
    Amber Trajectory     .crd
    Amber NetCDF         .nc
    Amber Restart        .rst7
    Amber NetCDF         .ncrst
    Charmm DCD           .dcd
    PDB                  .pdb
    Mol2                 .mol2
    Scripps              .binpos
    Gromacs              .trr
    SQM Input            .sqm
    ===================  =========

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> pt.write_traj("output/t.nc", traj, overwrite=True) # write to amber netcdf file

    >>> # write to multi pdb files (t.pdb.1, t.pdb.2, ...)
    >>> pt.write_traj("output/t.pdb", traj, overwrite=True, options='multi')

    >>> # write all frames to single pdb file and each frame is seperated by "MODEL" word
    >>> pt.write_traj("output/t.pdb", traj, overwrite=True, options='model')

    >>> # write to DCD file
    >>> pt.write_traj("output/test.dcd", traj, overwrite=True)

    >>> # write to netcdf file from 3D numpy array, need to provide Topology
    >>> xyz = traj.xyz
    >>> top = traj.top
    >>> pt.write_traj("output/test_xyz.nc", xyz, top=traj.top, overwrite=True)
    >>> pt.write_traj("output/test_xyz.nc", xyz, top=traj.top, overwrite=True)
    """
    from .frame import Frame
    from .trajs.Trajout import Trajout

    _top = _get_topology(traj, top)
    if _top is None:
        raise ValueError("must provide Topology")

    if traj is None or _top is None:
        raise ValueError("Need non-empty traj and top files")

    if not isinstance(traj, np.ndarray):
        with Trajout(filename=filename,
                     top=_top,
                     overwrite=overwrite,
                     options=options) as trajout:
            if isinstance(traj, Frame):
                if frame_indices is not None:
                    raise ValueError(
                        "frame indices does not work with single Frame")
                trajout.write(0, traj)
            else:
                if frame_indices is not None:
                    if isinstance(traj, (list, tuple, Frame)):
                        raise NotImplementedError(
                            "must be Trajectory or TrajectoryIterator instance")
                    for idx, frame in enumerate(traj.iterframe(
                        frame_indices=frame_indices)):
                        trajout.write(idx, frame)

                else:
                    for idx, frame in enumerate(iterframe_master(traj)):
                        trajout.write(idx, frame)
    else:
        # is ndarray, shape=(n_frames, n_atoms, 3)
        # create frame iterator
        xyz = np.asarray(traj)
        if not xyz.flags.c_contiguous:
            xyz = np.ascontiguousarray(xyz)
        _frame_indices = range(
            xyz.shape[0]) if frame_indices is None else frame_indices
        fi = iterframe_from_array(xyz, _top.n_atoms, _frame_indices)

        with Trajout(filename=filename,
                     top=_top,
                     overwrite=overwrite,
                     options=options) as trajout:

            for idx, frame in enumerate(fi):
                trajout.write(idx, frame)


def write_parm(filename=None, top=None, format='amberparm'):
    parm = ParmFile()
    parm.writeparm(filename=filename, top=top, format=format)


def load_topology(filename, more_options=''):
    """load Topology from a filename or from url or from ParmEd object. Adapted from cpptraj doc.

    Parameters
    ----------
    filename : str, Amber prmtop, pdb, mol2, psf, cif, gromacs topology, sdf, tinker formats
    more_options : cpptraj options.
        if filename is a pdb file, more_options = {'pqr', 'noconnect'}.
        pqr     : Read atomic charge/radius from occupancy/B-factor columns.
        noconect: Do not read CONECT records if present.

    Notes
    -----
    if cpptraj/pytraj does not support specific file format, you still can convert to PDB
    file. cpptraj will do the bond guess based on distance.

    Examples
    --------
    >>> import pytraj as pt
    >>> # from a filename
    >>> pt.load_topology("data/tz2.ortho.parm7")
    <Topology: 5293 atoms, 1704 residues, 1692 mols, PBC with box type = ortho>

    >>> # from url
    >>> pt.load_topology("http://ambermd.org/tutorials/advanced/tutorial1/files/polyAT.pdb")
    <Topology: 438 atoms, 438 residues, 438 mols, non-PBC>

    >>> # from ParmEd object
    >>> import parmed as pmd
    >>> parm = pmd.load_file('m2-c1_f3.mol2') # doctest: +SKIP
    >>> top = pt.load_topology(parm) # doctest: +SKIP

    >>> # read with more_options
    >>> pt.load_topology('1KX5.pdb', 'bondsearch 0.2') # doctest: +SKIP
    """
    top = Topology()

    # always read box info from pdb
    more_options = ' '.join(('readbox', more_options))

    if isinstance(filename, string_types):
        if filename.startswith('http://') or filename.startswith('https://'):
            import parmed as pmd
            return load_ParmEd(pmd.load_file(filename))
        else:
            parm = ParmFile()
            set_error_silent(True)
            parm.readparm(filename=filename,
                          top=top,
                          more_options=more_options)
            set_error_silent(False)
    else:
        raise ValueError('filename must be a string')

    if top.n_atoms == 0:
        raise RuntimeError(
            'n_atoms = 0: make sure to load correct filename '
            'or load supported topology (pdb, amber parm, psf, ...)')
    return top

# creat alias
read_parm = load_topology


def loadpdb_rcsb(pdbid):
    """load pdb file from rcsb website

    Parameters
    ----------
    pdbid : str

    Examples
    --------
        io.loadpdb_rcsb("2KOC") # popular RNA hairpin
    """
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid.upper()
    txt = urlopen(url).read()
    fname = "/tmp/tmppdb.pdb"
    with open(fname, 'w') as fh:
        if PY3:
            txt = txt.decode()
        fh.write(txt)
    traj = load(fname)
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


@ensure_exist
def load_single_frame(filename=None, top=None, index=0):
    """load a single Frame"""
    return iterload(filename, top)[index]


load_frame = load_single_frame

# creat alias
save = write_traj
save_traj = write_traj


def get_coordinates(iterable,
                    autoimage=None,
                    rmsfit=None,
                    mask=None,
                    frame_indices=None):
    '''return 3D-ndarray coordinates of `iterable`, shape=(n_frames, n_atoms, 3). This method is more memory
    efficient if use need to perform autoimage and rms fit to reference before loading all coordinates
    from disk.

    This method is good (fast, memory efficient) if you just want to get raw numpy array
    to feed to external package, such as sciki-learn, ...

    Parameters
    ----------
    iterable : could be anything that produces Frame when iterating 
               (Trajectory, TrajectoryIterator, FrameIterator, Python's generator, ...)

    Notes
    -----
    - if using both ``autoimage`` and ``rmsfit``, autoimage will be always processed before doing rmsfit.
    - You will get faster speed if ``iterable`` has attribute ``n_frames`` 

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2')
    >>> # load to out-of-core pytraj.TrajectoryIterator
    >>> traj = pt.iterload(fn, tn)
    >>> traj.n_frames, traj.n_atoms
    (10, 5293)

    >>> # simple
    >>> xyz = pt.get_coordinates(traj)  # same as traj.xyz
    >>> xyz.shape
    (10, 5293, 3)

    >>> # load coordinates to memory and doing autoimage
    >>> xyz = pt.get_coordinates(traj, autoimage=True)
    >>> xyz.shape
    (10, 5293, 3)

    >>> # load coordinates to memory for given mask
    >>> xyz = pt.get_coordinates(traj, mask='@CA')
    >>> xyz.shape
    (10, 12, 3)

    >>> # load coordinates of specific frame to memory and doing autoimage
    >>> xyz = pt.get_coordinates(traj, autoimage=True, frame_indices=[3, 6, 2, 5])
    >>> xyz.shape
    (4, 5293, 3)

    >>> # create frame iterator with some given cpptraj's commands
    >>> fi = pt.create_pipeline(traj, ['autoimage', 'rms', 'center :1-6 origin'])
    >>> xyz = pt.get_coordinates(fi)
    >>> xyz.shape
    (10, 5293, 3)

    >>> # make your own out-of-core method
    >>> def my_method(traj):
    ...     for frame in traj:
    ...         frame.xyz += 2.
    ...         yield frame
    >>> fi = my_method(traj)
    >>> fi.__class__.__name__
    'generator'
    >>> xyz = pt.get_coordinates(fi)
    >>> xyz.shape
    (10, 5293, 3)
    '''
    has_any_iter_options = any(
        x is not None for x in (autoimage, rmsfit, mask, frame_indices))
    # try to iterate to get coordinates
    if isinstance(iterable, (Trajectory, TrajectoryIterator)):
        fi = iterable.iterframe(autoimage=autoimage,
                                rmsfit=rmsfit,
                                mask=mask,
                                frame_indices=frame_indices)
    else:
        if has_any_iter_options:
            raise ValueError(
                'only support autoimage, rmsfit or mask for Trajectory and TrajectoryIterator')
        fi = iterframe_master(iterable)
    if hasattr(fi, 'n_frames') and hasattr(fi, 'n_atoms'):
        # faster
        n_frames = fi.n_frames
        shape = (n_frames, fi.n_atoms, 3)
        arr = np.empty(shape, dtype='f8')
        for idx, frame in enumerate(fi):
            # real calculation
            arr[idx] = frame.xyz
        return arr
    else:
        # slower
        return np.array([frame.xyz.copy()
                         for frame in iterframe_master(iterable)],
                        dtype='f8')
