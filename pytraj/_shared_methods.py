# distutils: language = c++
#
from pytraj.Frame import Frame
from pytraj.AtomMask import AtomMask
from pytraj.trajs.Trajout import Trajout
from pytraj.externals.six import string_types, set
from pytraj.utils import _import_numpy
from pytraj.utils.check_and_assert import is_frame_iter, is_chunk_iter
from pytraj.frameiter import FrameIter

__all__ = ['_savetraj', 'iterframe_master', '_xyz', 'my_str_method', '_tolist',
           '_box']


def _savetraj(self,
              filename="",
              format='unknown',
              overwrite=False, *args, **kwd):
    if format == 'unknown':
        # convert to "UNKNOWN_TRAJ"
        format = format.upper() + "_TRAJ"
    else:
        format = format.upper()

    with Trajout(filename=filename,
                 top=self.top,
                 format=format,
                 overwrite=overwrite, *args, **kwd) as trajout:
        for idx, frame in enumerate(self):
            trajout.write(idx, frame, self.top)



def _get_temperature_set(self):
    return set(self.temperatures)


def _xyz(self):
    """return a copy of xyz coordinates (wrapper of ndarray, shape=(n_frames, n_atoms, 3)
    We can not return a memoryview since Trajectory is a C++ vector of Frame object

    Notes
    -----
        read-only
    """
    import numpy as np
    n_frames = self.n_frames
    n_atoms = self.n_atoms

    myview = np.empty((n_frames, n_atoms, 3), dtype='f8')

    if self.n_atoms == 0:
        raise NotImplementedError("need to have non-empty Topology")
    for i, frame in enumerate(self):
        myview[i] = frame.buffer2d
    return myview


def _tolist(self):
    """return flatten list for traj-like object"""
    from itertools import chain
    return [frame.tolist() for frame in self]


def my_str_method(self):
    name = "pytraj." + self.__class__.__name__
    top_str = self.top.__str__()
    tmps = """<%s, %s frames, include:\n%s>
           """ % (
        name, self.n_frames, top_str, )
    return tmps


def _frame_iter(self, start=0, stop=-1, stride=1, mask=None):
    """iterately get Frames with start, stop, stride 
    Parameters
    ---------
    start : int (default = 0)
    stop : int (default = max_frames - 1)
    stride : int
    mask : str or array of interger
    """
    frame = Frame(self.n_atoms)

    if stop == -1:
        _end = self.n_frames
    else:
        _end = stop

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
                    raise 'TypeError'
            frame2 = Frame(frame, atm)
            yield frame2
        else:
            yield frame
        i += stride


def iterframe_master(obj):
    """try to return a frame iterator

    Parameters
    ----------
    obj : Trajectory or TrajectoryIterator or FrameIter or a list of frames, or a list of
    Trajectory

    Examples
    --------
    >>> from pytraj import iterframe_master
    >>> # create frame itrator from a TrajectoryIterator
    >>> for frame in iterframe_master(traj): pass

    >>> # create frame itrator from a list of trajs
    >>> for frame in iterframe_master([traj, traj]): pass

    >>> # create frame iterator from a list of traj and frame
    >>> for frame in iterframe_master([traj, traj[0]]): pass

    >>> # create frame iterator from a list of frame iterator and chunk iterator
    >>> for frame in iterframe_master([traj.iterframe(), traj.iterchunk()]): pass
    """

    is_frame_iter_but_not_master = (
        is_frame_iter(obj) and not obj.__name__ is 'iterframe_master')
    if isinstance(obj, Frame):
        yield obj
    elif isinstance(obj, FrameIter):
        for frame in obj:
            yield frame
    elif hasattr(obj, 'n_frames') or is_frame_iter_but_not_master:
        # traj-like or frame_iter or _frame_iter
        for frame in obj:
            yield frame
    else:
        try:
            # list, tuple, TrajinList, iterchunk
            for traj_obj in obj:
                if isinstance(traj_obj, Frame):
                    frame = traj_obj
                    yield frame
                elif is_chunk_iter(traj_obj):
                    for _traj in traj_obj:
                        for frame in _traj:
                            yield frame
                else:
                    for frame in traj_obj:
                        yield frame
        except:
            raise ValueError("can not convert to Frame")


def _box(self):
    import numpy as np
    boxarr = np.empty(self.n_frames * 6,
                      dtype=np.float64).reshape(self.n_frames, 6)

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
