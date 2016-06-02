from __future__ import absolute_import
import os
import numpy as np

from .externals.six import string_types, PY3
from .serialize import to_pickle, read_pickle, to_json, read_json
from .datafiles.load_samples import load_sample_data
from .datafiles.load_cpptraj_file import load_cpptraj_file
from .shared_methods import iterframe_master
from .cyutils import _fast_iterptr as iterframe_from_array
from .c_options import set_error_silent
from .get_common_objects import get_topology
from .topology import Topology, ParmFile
from .trajectory import Trajectory
from .trajectory_iterator import TrajectoryIterator

from .externals.load_other_packages import load_ParmEd

from .decorators import ensure_exist
from .core.c_core import _load_batch

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
           'write_parm',
           'save',
           'write_traj',
           'read_pickle',
           'read_json',
           'to_pickle',
           'to_json', ]


def load(filename, top=None, frame_indices=None, mask=None, stride=None):
    """try loading and returning appropriate values. See example below.

    Parameters
    ----------
    filename : str, Trajectory filename
    top : Topology filename or a Topology
    frame_indices : {None, array_like}, default None
        only load frames with given number given in frame_indices
    stride : {None, int}, default None
        if given, frame will be skip every `stride`.
        Note: if bot frame_indices and stride are given, `frame_indices` will be ignored.
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

    >>> # load mol2 file
    >>> traj = pt.load('traj.mol2') # # doctest: +SKIP

    >>> # load pdb file
    >>> traj = pt.load('traj.pdb') # doctest: +SKIP

    >>> # load given frame numbers
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2_dry')
    >>> traj = pt.load(fn, top=tn, frame_indices=[0, 3, 5, 12, 20])
    >>> traj = pt.load(fn, top=tn, frame_indices=[0, 3, 5, 12, 20], mask='!@H=')

    >>> # load with frame slice
    >>> traj = pt.load(fn, tn, frame_indices=slice(0, 10, 2))
    >>> # which is equal to:
    >>> traj = pt.load(fn, tn, frame_indices=range(0, 10, 2))

    >>> # load with stride
    >>> traj = pt.load(fn, tn)
    >>> traj.n_frames
    101
    >>> traj = pt.load(fn, tn, stride=5)
    >>> traj.n_frames
    21

    >>> # load with stride for more than one filename
    >>> traj = pt.load([fn, fn], tn, stride=5)
    >>> traj.n_frames
    42
    >>> traj.n_atoms
    223

    >>> # load with stride for more than one filename, and with mask
    >>> traj = pt.load([fn, fn], tn, stride=5, mask='@CA')
    >>> traj.n_frames
    42
    >>> traj.n_atoms
    12
    """
    # load to TrajectoryIterator object first
    # do not use frame_indices_ here so we can optimize the slicing speed
    traj = load_traj(filename, top, stride=stride)

    # do the slicing and other things if needed.
    if stride is not None:
        if mask is None:
            return traj[:]
        else:
            return traj[mask]
    else:
        frame_indices_ = frame_indices
        if isinstance(frame_indices_, tuple):
            frame_indices_ = list(frame_indices_)
        if frame_indices_ is None and mask is None:
            # load all
            return traj[:]
        elif frame_indices_ is None and mask is not None:
            # load all frames with given mask
            # eg: traj['@CA']
            return traj[mask]
        elif frame_indices_ is not None and mask is None:
            # eg. traj[[0, 3, 7]]
            return traj[frame_indices_]
        else:
            # eg. traj[[0, 3, 7], '@CA']
            return traj[frame_indices_, mask]


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
    stride : {None, int}, default None
        if not None, trajectories will be strided.
        Note: if both stride and frame_slice are not None, frame_slice will be ignored

    Returns
    -------
    pytraj.TrajectoryIterator

    Examples
    --------
    >>> import pytraj as pt
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
    >>> traj = pt.iterload(['traj0.nc', 'traj1.nc'], '2koc.parm7', frame_slice=[(0, 10), (0, -1)]) # doctest: +SKIP

    >>> # use stride, skip every 2 frames
    >>> from pytraj.testing import get_remd_fn
    >>> filenames, topology_filename = get_remd_fn('remd_ala2')
    >>> [fn.split('/')[-1] for fn in filenames]
    ['rem.nc.000', 'rem.nc.001', 'rem.nc.002', 'rem.nc.003']
    >>> traj = pt.iterload(filenames, topology_filename, stride=2)

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

    if 'stride' in kwd:
        ts._load(filename, stride=kwd['stride'])
    elif 'frame_slice' in kwd:
        ts._load(filename, frame_slice=kwd['frame_slice'])
    else:
        ts._load(filename)

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
    from pytraj.core.c_core import CpptrajState, Command

    state = CpptrajState()

    # add keyword 'remdtraj' to trick cpptraj
    trajin = ' '.join(('trajin', filename, 'remdtraj remdtrajtemp', str(T)))
    if isinstance(top, string_types):
        top = load_topology(top)
    else:
        top = top
    state.data.add('topology', 'remdtop')
    # set topology
    state.data['remdtop']._top = top

    # load trajin
    with Command() as cm:
        cm.dispatch(state, trajin)
        cm.dispatch(state, 'loadtraj name remdtraj')

    # state.data.remove_set(state.data['remdtop'])
    traj = state.data[-1]

    # assign state to traj to avoid memory free
    traj._base = state
    return traj


def load_remd(filename, top=None, T="300.0"):
    return iterload_remd(filename, top, T)[:]


def write_traj(filename="",
               traj=None,
               format='infer',
               top=None,
               frame_indices=None,
               overwrite=False,
               crdinfo=None,
               options=""):
    """write Trajectory-like or iterable object to trajectory file

    Parameters
    ----------
    filename : str
    traj : Trajectory-like or iterator that produces Frame or 3D ndarray with shape=(n_frames, n_atoms, 3)
    format : str, default 'infer'
        if 'inter', detect format based on extension. If can not detect, use amber mdcdf format.
    top : Topology, optional, default: None
    frame_indices: array-like or iterator that produces integer, default: None
        If not None, only write output for given frame indices
    overwrite: bool, default: False
    crdinfo : None or dict, default None
        if None, try to get info from traj._crdinfo (if traj has _crdinfo)
        if given, use it. `crdinfo` needed to pass if you want to write force/velocity (netcdf)
    options : str, additional cpptraj keywords

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

    'options' for writing to pdb format (cptraj manual)::

        dumpq:       Write atom charge/GB radius in occupancy/B-factor columns (PQR format)."
        parse:       Write atom charge/PARSE radius in occupancy/B-factor columns (PQR format)."
        vdw:         Write atom charge/VDW radius in occupancy/B-factor columns (PQR format)."
        pdbres:      Use PDB V3 residue names."
        pdbatom:     Use PDB V3 atom names."
        pdbv3:       Use PDB V3 residue/atom names."
        teradvance:  Increment record (atom) # for TER records (default no)."
        terbyres:    Print TER cards based on residue sequence instead of molecules."
        model:       Write to single file separated by MODEL records."
        multi:       Write each frame to separate files."
        chainid <c>: Write character 'c' in chain ID column."
        sg <group>:  Space group for CRYST1 record, only if box coordinates written."
        include_ep:  Include extra points."
        conect:      Write CONECT records using bond information.");

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

    >>> # you can make a fake Topology to write xyz coordinates too
    >>> n_atoms = xyz.shape[1]
    >>> top2 = pt.tools.make_fake_topology(n_atoms)
    >>> pt.write_traj("output/test_xyz_fake_top.nc", xyz, top=top2, overwrite=True)

    'options' for writing to amber netcdf format (cptraj manual)::

        remdtraj: Write temperature to trajectory (makes REMD trajectory)."
        velocity: Write velocities to trajectory."
        force: Write forces to trajectory.");

    'options' for writing to amber netcdf restart format(cptraj manual)::

        novelocity: Do not write velocities to restart file."
        notime:     Do not write time to restart file."
        remdtraj:   Write temperature to restart file."
        time0:      Time for first frame (default 1.0)."
        dt:         Time step for subsequent frames, t=(time0+frame)*dt; (default 1.0)");
        keepext     Keep filename extension; write '<name>.<num>.<ext>' instead (example: myfile.1.rst7)

    'options' for writing to mol2 format (cptraj manual)::

        single   : Write to a single file."
        multi    : Write each frame to a separate file."
        sybyltype: Convert Amber atom types (if present) to SYBYL types.");

    'options'  for other formats::

        please check http://ambermd.org/doc12/Amber15.pdf
    """
    from .frame import Frame
    from .c_traj.c_trajout import TrajectoryWriter

    _top = get_topology(traj, top)
    if _top is None:
        raise ValueError("must provide Topology")

    if traj is None or _top is None:
        raise ValueError("Need non-empty traj and top files")

    if not isinstance(traj, np.ndarray):
        with TrajectoryWriter(filename=filename,
                     top=_top,
                     format=format,
                     overwrite=overwrite,
                     options=options) as trajout:
            if isinstance(traj, Frame):
                if frame_indices is not None:
                    raise ValueError(
                        "frame indices does not work with single Frame")
                trajout.write(traj)
            else:
                if frame_indices is not None:
                    if isinstance(traj, (list, tuple, Frame)):
                        raise NotImplementedError(
                            "must be Trajectory or TrajectoryIterator instance")
                    for frame in traj.iterframe(frame_indices=frame_indices):
                        trajout.write(frame)

                else:
                    for frame in iterframe_master(traj):
                        trajout.write(frame)
    else:
        # is ndarray, shape=(n_frames, n_atoms, 3)
        # create frame iterator
        xyz = np.asarray(traj)
        if not xyz.flags.c_contiguous:
            xyz = np.ascontiguousarray(xyz)
        _frame_indices = range(xyz.shape[
            0]) if frame_indices is None else frame_indices
        fi = iterframe_from_array(xyz, _top.n_atoms, _frame_indices, _top)

        if crdinfo is None:
            if hasattr(traj, '_crdinfo'):
                crdinfo = traj._crdinfo
            else:
                crdinfo = dict()
        else:
            crdinfo = crdinfo

        with TrajectoryWriter(filename=filename,
                     top=_top,
                     crdinfo=crdinfo,
                     overwrite=overwrite,
                     options=options) as trajout:

            for frame in fi:
                trajout.write(frame)


def write_parm(filename=None, top=None, format='amberparm', overwrite=False):
    if os.path.exists(filename) and not overwrite:
        raise RuntimeError('{0} exists, must set overwrite=True'.format(
            filename))
    parm = ParmFile()
    parm.writeparm(filename=filename, top=top, format=format)


def load_topology(filename, option=''):
    """load Topology from a filename or from url or from ParmEd object. Adapted from cpptraj doc.

    Parameters
    ----------
    filename : str, Amber prmtop, pdb, mol2, psf, cif, gromacs topology, sdf, tinker formats
    option : cpptraj options.
        if filename is a pdb file, option = {'pqr', 'noconnect'}.

        - pqr     : Read atomic charge/radius from occupancy/B-factor columns.
        - noconect: Do not read CONECT records if present.

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

    >>> # read with option
    >>> pt.load_topology('1KX5.pdb', 'bondsearch 0.2') # doctest: +SKIP
    """
    top = Topology()

    # always read box info from pdb
    option = ' '.join(('readbox', option))

    if isinstance(filename, string_types):
        parm = ParmFile()
        set_error_silent(True)
        parm.readparm(filename=filename,
                      top=top,
                      option=option)
        set_error_silent(False)
    else:
        raise ValueError('filename must be a string')

    if top.n_atoms == 0:
        raise RuntimeError(
            'n_atoms = 0: make sure to load correct Topology filename '
            'or load supported topology (pdb, amber parm, psf, ...)')
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
    url = 'http://files.rcsb.org/download/{}.pdb'.format(pdbid.upper())
    return _make_traj_from_remote_file(url)

def load_url(url):
    """
    
    versionadded: 1.0.7
    """
    return _make_traj_from_remote_file(url)

def _make_traj_from_remote_file(remote_file):
    import tempfile

    fd, fname = tempfile.mkstemp()
    txt = urlopen(remote_file).read()
    with open(fname, 'w') as fh:
        if PY3:
            txt = txt.decode()
        fh.write(txt)

    return load(fname)

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
save_traj = write_traj


def save(filename, obj, *args, **kwd):
    '''an universal method

    Parameters
    ----------
    filename : output filename
    obj : Topology or Trajetory-like
        if Topology, write a new Topology to disk
        if Trajetory-like, write a trajectory to disk
    *args, **kwd: additional options

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_ala3()
    >>> isinstance(traj, pt.TrajectoryIterator)
    True
    >>> top = traj.top
    >>> isinstance(top, pt.Topology)
    True
    >>> # save Topology to a new Topology
    >>> pt.save('output/prmtop', top, overwrite=True)
    >>> isinstance(pt.load_topology('output/prmtop'), pt.Topology)
    True
    >>> # save TrajectoryIterator to disk
    >>> pt.save('output/traj.nc', traj, overwrite=True)
    >>> isinstance(pt.iterload('output/traj.nc', traj.top), pt.TrajectoryIterator)
    True

    See also
    --------
    write_traj
    write_parm
    '''
    if isinstance(obj, Topology):
        write_parm(filename, obj, *args, **kwd)
    else:
        write_traj(filename, obj, *args, **kwd)


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
    >>> fi = pt.pipe(traj, ['autoimage', 'rms', 'center :1-6 origin'])
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
        return np.array(
            [frame.xyz.copy() for frame in iterframe_master(iterable)],
            dtype='f8')

def load_batch(traj, txt):
    '''perform calculation for traj with cpptraj's batch style. This is for internal use.

    Parameters
    ----------
    traj : pytraj.TrajectoryIterator
    txt : text or a list of test
        cpptraj's commands

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> text = """
    ... autoimage
    ... radgyr @CA nomax
    ... molsurf !@H=
    ... """
    >>> state = pt.load_batch(traj, text)
    >>> state = state.run()
    >>> state.data
    <pytraj.datasets.CpptrajDatasetList - 3 datasets>

    >>> # raise if not TrajectoryIterator
    >>> traj2 = pt.Trajectory(xyz=traj.xyz, top=traj.top)
    >>> not isinstance(traj2, pt.TrajectoryIterator)
    True
    >>> pt.load_batch(traj2, text)
    Traceback (most recent call last):
        ...
    ValueError: only support TrajectoryIterator
    '''
    if not isinstance(traj, TrajectoryIterator):
        raise ValueError('only support TrajectoryIterator')
    return _load_batch(txt, traj=traj)
