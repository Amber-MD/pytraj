# distutils: language = c++
#
cimport cython
from pytraj.Frame cimport _Frame, Frame
from pytraj.Trajectory cimport Trajectory
from pytraj.AtomMask cimport AtomMask
from pytraj.trajs.Trajout import Trajout
from pytraj.externals.six import string_types
from pytraj.compat import set
from pytraj.utils import _import_numpy
from pytraj.exceptions import PytrajMemviewError, PytrajConvertError
from pytraj.utils.check_and_assert import is_frame_iter, is_chunk_iter
from pytraj._xyz import XYZ

__all__ = ['_savetraj', '_frame_iter_master', '_xyz', 'my_str_method',
           '_tolist', '_box_to_ndarray']

def _savetraj(self, filename="", fmt='unknown', overwrite=False, *args, **kwd):
    if fmt == 'unknown':
        # convert to "UNKNOWN_TRAJ"
        fmt = fmt.upper() + "_TRAJ"
    else:
        fmt = fmt.upper()

    with Trajout(filename=filename, top=self.top, fmt=fmt, 
                 overwrite=overwrite, *args, **kwd) as trajout:
        for idx, frame in enumerate(self):
            trajout.writeframe(idx, frame, self.top)

def _split_and_write_traj(self, n_chunks=None, root_name="trajx", *args, **kwd):
    chunksize = self.n_frames // n_chunks
    for idx, traj in enumerate(self.chunk_iter(chunksize=chunksize)):
        fname = ".".join((root_name, str(idx), 'nc'))
        traj.save(fname, *args, **kwd)

def _get_temperature_set(self):
    return set(self.temperatures) 

def _xyz(self):
    """return a copy of xyz coordinates (wrapper of ndarray, shape=(n_frames, n_atoms, 3)
    We can not return a memoryview since Trajectory is a C++ vector of Frame object

    Notes
    -----
        read-only
    """
    cdef bint has_numpy
    cdef int i
    cdef int n_frames = self.n_frames
    cdef int n_atoms = self.n_atoms
    cdef Frame frame

    has_numpy, np = _import_numpy()
    myview = np.empty((n_frames, n_atoms, 3), dtype='f8')

    if self.n_atoms == 0:
        raise NotImplementedError("need to have non-empty Topology")
    if has_numpy:
        for i, frame in enumerate(self):
            myview[i] = frame.buffer2d
        return XYZ(myview)
    else:
        raise NotImplementedError("must have numpy")

def _tolist(self):
    """return flatten list for traj-like object"""
    from itertools import chain
    return [frame.tolist() for frame in self]

def my_str_method(self):
    name = "pytraj." + self.__class__.__name__
    top_str = self.top.__str__()
    tmps = """<%s with %s frames: %s>
           """ % (
            name, self.size, top_str,
            )
    return tmps

def _frame_iter(self, int start=0, int stop=-1, int stride=1, mask=None):
    """iterately get Frames with start, stop, stride 
    Parameters
    ---------
    start : int (default = 0)
    stop : int (default = max_frames - 1)
    stride : int
    mask : str or array of interger
    """
    cdef int i
    cdef Frame frame = Frame(self.n_atoms)
    cdef Frame frame2
    cdef AtomMask atm
    cdef int _end
    cdef int[:] int_view

    if stop == -1:
        _end = <int> self.n_frames
    else:
        _end = stop + 1

    i = start
    while i < _end:
        frame = self[i]
        if mask is not None:
            if isinstance(mask, string_types):
                atm = self.top(mask)
            else:
                try:
                    atm = AtomMask()
                    atm.add_selected_indices(mask)
                except TypeError:
                    raise PytrajMemviewError()
            frame2 = Frame(atm.n_atoms)
            frame2.thisptr.SetCoordinates(frame.thisptr[0], atm.thisptr[0])
            yield frame2
        else:
            yield frame
        i += stride

def _frame_iter_master(obj):
    """try to return frame iterator object

    obj : could be Trajectory, TrajectoryIterator, TrajinList, frame_iter object
          could be a list of trajs, ...

    Examples
    --------
    >>> from pytraj import io
    >>> from pytraj import Frame
    >>> traj = io.load_sample_data('tz2')
    >>> from pytraj._shared_methods import _frame_iter_master
    >>> for frame in _frame_iter_master(traj): assert isinstance(frame, Frame) == True
    >>> for frame in _frame_iter_master([traj,]): assert isinstance(frame, Frame) == True
    >>> for frame in _frame_iter_master([traj, traj[0]]): assert isinstance(frame, Frame) == True
    >>> for frame in _frame_iter_master([traj.frame_iter(), traj.chunk_iter()]): 
    >>>     assert isinstance(frame, Frame) == True
    """
    cdef Frame frame
    cdef object traj_obj
    cdef Trajectory _traj

    is_frame_iter_but_not_master = (is_frame_iter(obj) and not obj.__name__ is '_frame_iter_master')
    if isinstance(obj, Frame):
        yield obj
    elif hasattr(obj, 'n_frames') or is_frame_iter_but_not_master:
        # traj-like or frame_iter or _frame_iter
        for frame in obj:
            yield frame
    else:
        try:
            # list, tuple, TrajinList, chunk_iter
            for traj_obj in obj:
                if isinstance(traj_obj, Frame):
                    frame = <Frame> traj_obj
                    yield frame
                elif is_chunk_iter(traj_obj):
                    for _traj in traj_obj:
                        for frame in _traj:
                            yield frame
                else:
                    for frame in traj_obj:
                        yield frame
        except:
            raise PytrajConvertError("can not convert to Frame")

def _box_to_ndarray(self): 
    cdef Frame frame
    cdef int i

    _, np = _import_numpy()
    boxarr = np.empty(self.n_frames * 6, dtype=np.float64).reshape(self.n_frames, 6)

    # Note: tried `enumerate` but got wrong result.
    # --> use old fashion
    i = 0
    for frame in self:
        boxarr[i] = frame.box.to_ndarray()
        i += 1
    return boxarr

if __name__ == '__main__':
    import doctest
    doctest.testmod()
