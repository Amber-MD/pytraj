from pytraj.externals.six import string_types, PY3
from pytraj.base import *
from pytraj.iterload import _iterload
from pytraj.TrajReadOnly import TrajReadOnly
from pytraj.FrameArray import FrameArray
from pytraj.Topology import Topology
from pytraj.utils.check_and_assert import make_sure_exist
try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

__all__ = ['load', 'iterload', 'writetraj', 'readparm', 'writeparm']

def load(*args, **kwd):
    """try loading and returning appropriate values"""

    if 'filename' in kwd.keys():
        make_sure_exist(kwd['filename'])
    else:
        make_sure_exist(args[0])

    if len(args) + len(kwd) == 1:
        # loading only Topology
        return readparm(*args, **kwd)
    else:
        # load traj
        return loadtraj(*args, **kwd)

def iterload(filename=None, top=Topology(), start=0, chunk=None):
    """iterload(filename=None, top=Topology(), start=0, chunk=None)"""
    if not isinstance(top, Topology):
        # string
        top = Topology(top)
    #filename = filename.encode("UTF-8")
    return _iterload(top, filename, start, chunk)

def loadtraj(filename=None, top=Topology(), readonly=True, indices=None):
    """load(filename=None, top=Topology(), readonly=True)"""
    #filename = filename.encode("UTF-8")
    if not isinstance(top, Topology):
        # string
        #top = top.encode("UTF-8")
        top = Topology(top)
    if readonly:
        ts = TrajReadOnly()
    else:
        ts = FrameArray()
    # TODO : use indices
    if indices is not None:
        raise RuntimeError("not yet supported")
    ts.load(filename, top)
    return ts

def writetraj(filename="", traj=None, top=None, 
              fmt='UNKNOWN_TRAJ', indices=None,
              overwrite=False):
    """writetraj(filename="", traj=None, top=None, 
              ftm='UNKNOWN_TRAJ', indices=None):
    """
    # TODO : support list (tuple) of FrameArray, TrajReadOnly or 
    # list of filenames
    #filename = filename.encode("UTF-8")

    if fmt == 'unknown':
        fmt = fmt.upper() + "_TRAJ"
    else:
        fmt = fmt.upper()

    if not traj or not top:
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
                for idx, frame in enumerate(traj2):
                    trajout.writeframe(idx, frame, top)
            else:
                if isinstance(traj2, (list, tuple)):
                    raise NotImplementedError("must be FrameArray or TrajReadOnly instance")
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
    return Topology(filename)

def loadpdb_rcsb(pdbid):
    # TODO : use tempfile
    """load pdb file from rcsb website
    Return : Topology and FrameArray instance
    """

    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
    txt = urlopen(url).read()
    print ("test saving pdb")
    print (type(txt))
    with open("/tmp/._tmp", 'w') as fh:
        if PY3:
            fh.write(txt.decode())
        else:
            fh.write(txt)
    top = readparm("/tmp/._tmp")
    frames = load("/tmp/._tmp", top)
    return top, frames

def load_single_frame(frame=None, top=None):
    """load single Frame"""
    return load(frame, top)[0]
