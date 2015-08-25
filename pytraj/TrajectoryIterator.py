"""This is a thin wrapper of Trajin_Single
We need to sub-class Trajin_Single to use Trajectory
(we called Trajin_Single from Trajectory, so we can not call Trajectory back from 
Trajin_Single)
"""
from __future__ import absolute_import
import warnings
import os
import numpy as np

from .trajs.TrajectoryCpptraj import TrajectoryCpptraj
from ._action_in_traj import ActionTrajectory
from .compat import string_types, range
from ._shared_methods import _tolist, _split_and_write_traj
from .Topology import Topology
from .utils import is_int
from ._cyutils import get_positive_idx


__all__ = ['TrajectoryIterator', 'split_iterators']


def split_iterators(traj, n_chunks=1, start=0, stop=-1, stride=1, mask=None,
                    autoimage=False, rmsfit=None):
    return traj.split_iterators(n_chunks, start, stop, stride, mask, autoimage, rmsfit)


def _make_frame_slices(n_files, original_frame_slice):
    if isinstance(original_frame_slice, tuple):
        return [original_frame_slice for i in range(n_files)]
    elif isinstance(original_frame_slice, list):
        fs_len = len(original_frame_slice)
        if fs_len < n_files:
            new_list = original_frame_slice[:] + [(0, -1, 1)
                                                  for _ in range(fs_len, n_files)]
        elif fs_len == n_files:
            new_list = original_frame_slice
        else:
            raise ValueError(
                "len of frame_slice tuple-list must be smaller or equal number of files")
        return new_list
    else:
        raise ValueError(
            "must be a tuple of integer values or a list of tuple of integer values")


def _turn_to_list_with_rank(func):
    def inner(*args, **kwd):
        if 'rank' not in kwd:
            return list(func(*args, **kwd))
        else:
            return list(func(*args, **kwd))[kwd['rank']]
    inner.__doc__ = func.__doc__
    return inner


class TrajectoryIterator(TrajectoryCpptraj):

    def __init__(self, filename=None, top=None, *args, **kwd):
        self._force_load = False
        # use self.chunk to store `chunk` in chunk_iter
        # to deallocate memory
        self.chunk = None
        # same as self.chunk but for frame_iter
        self.frame = None
        # only allow to load <= 1000 Mb
        self._size_limit_in_MB = 1000
        super(TrajectoryIterator, self).__init__()

        if not top:
            self.top = Topology()
        elif isinstance(top, string_types):
            self.top = Topology(top)
        elif isinstance(top, Topology):
            self.top = top.copy()
        else:
            raise ValueError("Topology must be None/string/Topology")

        self.frame_slice_list = []

        if filename:
            if self.top.is_empty():
                raise ValueError('First argument is always a trajectory filename'
                                 ' or a list of filenames'
                                 'must have a non-empty Topology')
            self.load(filename, self.top, *args, **kwd)
        if not top and (args or kwd):
            warnings.warn('creating an empty TrajectoryIterator since does not '
                          'have Topology information. Ignore other arguments')

    def __iter__(self):
        for frame in super(TrajectoryIterator, self).__iter__():
            # we need to use `copy` to create different frame pointer
            # so [frame for frame in traj] will return a list of different ones
            yield frame.copy()

    def copy(self):
        """Very simple"""
        other = self.__class__()
        other.top = self.top.copy()
        other.load(self.filelist)

        return other

    def load(self, filename=None, top=None, frame_slice=(0, -1, 1)):
        """load trajectory/trajectories from filename/filenames
        with a single frame_slice or a list of frame_slice
        """
        if not top:
            _top = self.top
        else:
            _top = top

        if isinstance(filename, string_types) and os.path.exists(filename):
            super(TrajectoryIterator, self).load(filename, _top, frame_slice)
            self.frame_slice_list.append(frame_slice)
        elif isinstance(filename, (list, tuple)):
            filename_list = filename
            full_frame_slice = _make_frame_slices(
                len(filename_list), frame_slice)

            for fname, fslice in zip(filename_list, full_frame_slice):
                self.frame_slice_list.append(frame_slice)
                super(TrajectoryIterator, self).load(
                    fname, _top, frame_slice=fslice)
        elif isinstance(filename, string_types) and not os.path.exists(filename):
            from glob import glob
            flist = sorted(glob(filename))
            if not flist:
               raise ValueError("must provie a filename or list of filenames or file pattern")
            self.load(flist, top=top, frame_slice=frame_slice)
        else:
            raise ValueError("")

    @property
    def topology(self):
        """traditional name for Topology file"""
        return self.top

    @topology.setter
    def topology(self, newtop):
        self.top = newtop

    @property
    def coordinates(self):
        """return 3D numpy.ndarray, same as `TrajectoryIterator.xyz`
        """
        return self.xyz

    @property
    def _estimated_MB(self):
        """esimated MB of data will be loaded to memory
        """
        return self.n_frames * self.n_atoms * 3 * 8 / 1E6

    @property
    def xyz(self):
        size_in_MB = self._estimated_MB
        # check if larger than size_limit_in_MB
        if size_in_MB > self._size_limit_in_MB and not self._force_load:
            raise MemoryError("you are loading %s Mb, larger than size_limit %s Mb. "
                              "Please increase self._size_limit_in_MB or set self._force_load=True"
                              % (size_in_MB, self._size_limit_in_MB))
        return super(TrajectoryIterator, self).xyz

    def iterator_slice(self, start=0, stop=None, stride=None):
        """iterator slice"""
        from itertools import islice
        if stop is None:
            stop = self.n_frames
        return islice(self, start, stop, stride)

    def make_independent_iterators(self, n_iters):
        from itertools import tee
        return tee(self, n_iters)

    def frame_iter(self, start=0, stop=None, stride=1, mask=None,
                   autoimage=False, rmsfit=None):

        from .core.frameiter import FrameIter
        if mask is None:
            _top = self.top
        else:
            _top = self.top._get_new_from_mask(mask)

        if rmsfit is not None:
            if isinstance(rmsfit, tuple):
                assert len(rmsfit) <= 2, ("rmsfit must be a tuple of one (frame,) "
                                          "or two elements (frame, mask)")
                if len(rmsfit) == 1:
                    rmsfit = (rmsfit, '*')
            elif isinstance(rmsfit, int):
                rmsfit = (rmsfit, '*')
            else:
                raise ValueError("rmsfit must be a tuple or an integer")

            if is_int(rmsfit[0]):
                index = rmsfit[0]
                rmsfit = ([self[index], rmsfit[1]])

        # check how many frames will be calculated
        if stop is None or stop >= self.n_frames:
            stop = self.n_frames
        elif stop < 0:
            stop = get_positive_idx(stop, self.n_frames)

        n_frames = len(range(start, stop, stride))

        frame_iter_super = super(
            TrajectoryIterator, self).frame_iter(start, stop, stride)

        return FrameIter(frame_iter_super,
                         original_top=self.top,
                         new_top=_top,
                         start=start,
                         stop=stop,
                         stride=stride,
                         mask=mask,
                         autoimage=autoimage,
                         rmsfit=rmsfit,
                         is_trajiter=True,
                         n_frames=n_frames,
                         )

    def iterframe(self, *args, **kwd):
        return self.frame_iter(*args, **kwd)

    def iterchunk(self, *args, **kwd):
        return self.chunk_iter(*args, **kwd)

    def chunk_iter(self, chunksize=2, start=0, stop=-1,
                   autoimage=False,
                   rmsfit=None,
                   copy_top=False):
        """
        Parameters
        ----------
        chunk : int, default=2
            size of each chunk. Notes: final chunk's size might be changed
        start : int, default=0 (first frame)
        stop : int, default=-1 (last frame)
        autoimage : bool, default=False
        rmsfit : None | tuple/list of (reference frame, mask)

        Examples
        --------
            for chunk in trajiter.chunk_iter(100, autoimage=True, rmsfit=(ref0, '@CA'))
        """
        chunk = chunksize

        if rmsfit is not None:
            ref, mask_for_rmsfit = rmsfit
            need_align = True
        else:
            need_align = False
            ref, mask_for_rmsfit = None, None

        for chunk in super(TrajectoryIterator, self).chunk_iter(chunk,
                                                                start, stop, copy_top):
            # always perform autoimage before doing fitting
            # chunk is `Trajectory` object, having very fast `autoimage` and
            # `rmsfit` methods
            if autoimage:
                chunk.autoimage()
            if need_align:
                chunk.rmsfit(ref, mask_for_rmsfit)
            # free memory
            # if not, memory will be quicly accumulated.
            if self.chunk:
                self.chunk.__del__()

            self.chunk = chunk
            yield self.chunk

    @property
    def filename(self):
        '''return 1st filename in filelist. For testing only
        '''
        return self.filelist[0]

    @property
    def shape(self):
        return (self.n_frames, self.n_atoms, 3)

    def is_empty(self):
        return self.n_frames == 0

    @property
    def max_frames(self):
        return self.n_frames

    def tolist(self):
        return _tolist(self)

    def to_mutable_trajectory(self):
        return self[:]

    def split_and_write_traj(self, *args, **kwd):
        _split_and_write_traj(self, *args, **kwd)

    @_turn_to_list_with_rank
    def split_iterators(self, n_chunks=1, start=0, stop=-1, stride=1,
                        mask=None,
                        autoimage=False, rmsfit=None, **kwd):
        """simple splitting `self` to n_chunks FrameIter objects        

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> list(traj.split_iterators(n_chunks=4, mask='@CA'))
        [<pytraj.core.frameiter.FrameIter with start=0, stop=2, stride=1
         autoimage=False, rmsfit=None> ,
         <pytraj.core.frameiter.FrameIter with start=2, stop=4, stride=1
         autoimage=False, rmsfit=None> ,
         <pytraj.core.frameiter.FrameIter with start=4, stop=6, stride=1
         autoimage=False, rmsfit=None> ,
         <pytraj.core.frameiter.FrameIter with start=6, stop=10, stride=1
         autoimage=False, rmsfit=None> ]
        """
        from pytraj.tools import split_range

        assert 0 <= start <= self.n_frames, "0 <= start <= self.n_frames"

        if stop <= 0 or stop > self.n_frames:
            stop = self.n_frames

        if n_chunks == 1:
            yield self(start=start, stop=stop, stride=stride, mask=mask,
                       autoimage=autoimage, rmsfit=rmsfit)
        else:
            for (_start, _stop) in split_range(n_chunks=n_chunks,
                                               start=start, stop=stop):
                yield self.frame_iter(start=_start, stop=_stop, stride=stride,
                                      mask=mask,
                                      autoimage=autoimage,
                                      rmsfit=rmsfit)

    def to_numpy_traj(self):
        from pytraj import api
        return api.Trajectory(self)

    @property
    def temperatures(self):
        return np.array([frame.temperature for frame in self])
