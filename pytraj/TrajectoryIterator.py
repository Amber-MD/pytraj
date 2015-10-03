"""out-of-core TrajectoryIterator
"""
from __future__ import absolute_import
import warnings
import os
from glob import glob
import numpy as np
from .trajs.TrajectoryCpptraj import TrajectoryCpptraj
from .externals.six import string_types
from .externals.six.moves import range
from .Topology import Topology
from .Frame import Frame
from .utils import is_int
from ._cyutils import get_positive_idx
from .frameiter import FrameIter
from ._get_common_objects import _load_Topology
from .utils import split_range

__all__ = ['TrajectoryIterator', 'split_iterators']


def split_iterators(traj,
                    n_chunks=1,
                    start=0,
                    stop=-1,
                    stride=1,
                    mask=None,
                    autoimage=False,
                    rmsfit=None):
    return traj._split_iterators(n_chunks, start, stop, stride, mask, autoimage,
                                rmsfit, rank=-1)


def _make_frame_slices(n_files, original_frame_slice):
    if isinstance(original_frame_slice, tuple):
        return [original_frame_slice for i in range(n_files)]
    elif isinstance(original_frame_slice, list):
        fs_len = len(original_frame_slice)
        if fs_len < n_files:
            new_list = original_frame_slice[:] + [(
                0, -1, 1) for _ in range(fs_len, n_files)]
        elif fs_len == n_files:
            new_list = original_frame_slice
        else:
            raise ValueError(
                "len of frame_slice tuple-list must be smaller or equal number of files")
        return new_list
    else:
        raise ValueError(
            "must be a tuple of integer values or a list of tuple of integer values")


class TrajectoryIterator(TrajectoryCpptraj):
    def __init__(self, filename=None, top=None, *args, **kwd):
        '''out-of-core trajectory holder.

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.TrajectoryIterator('./traj.nc', 'input.parm7')

        >>> # user should always use :method:`pytraj.iterload` to load TrajectoryIterator
        >>> traj = pt.iterload(['remd.x.000', 'remd.x.001'], 'input.parm7')
        >>> # load another trajectory
        >>> traj.load('./remd.x.003')
        '''
        self._force_load = False
        # use self._chunk to store `chunk` in iterchunk
        # to deallocate memory
        self._chunk = None
        # only allow to load <= 1 GB
        self._size_limit_in_GB = 1
        super(TrajectoryIterator, self).__init__()

        if not top:
            self.top = Topology()
        elif isinstance(top, string_types):
            self.top = _load_Topology(top)
        elif isinstance(top, Topology):
            self.top = top.copy()
        else:
            raise ValueError("Topology must be None/string/Topology")

        self.frame_slice_list = []

        if filename:
            if self.top.is_empty():
                raise ValueError(
                    'First argument is always a trajectory filename'
                    ' or a list of filenames'
                    'must have a non-empty Topology')
            self.load(filename, self.top, *args, **kwd)
        if not top and (args or kwd):
            warnings.warn(
                'creating an empty TrajectoryIterator since does not '
                'have Topology information. Ignore other arguments')

        self.__dict__.update({
            'top': self.top,
            'top_filename': self.top.filename,
            'filelist': self.filelist
        })

    def __setstate__(self, state):
        self.__dict__ = state.copy()
        self.top = _load_Topology(state['top_filename'])
        self.load(state['filelist'], frame_slice=state['frame_slice_list'])

    def __getstate__(self):
        return self.__dict__

    def __iter__(self):
        '''do not make a frame copy here
        '''
        for frame in super(TrajectoryIterator, self).__iter__():
            yield frame

    def copy(self):
        '''return a deep copy
        '''
        other = self.__class__()
        other.top = self.top.copy()

        for fname, frame_slice in zip(self.filelist, self.frame_slice_list):
            other.load(fname, frame_slice=frame_slice)
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
                    fname, _top,
                    frame_slice=fslice)
        elif isinstance(filename,
                        string_types) and not os.path.exists(filename):
            flist = sorted(glob(filename))
            if not flist:
                raise ValueError(
                    "must provie a filename or list of filenames or file pattern")
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
    def _estimated_GB(self):
        """esimated GB of data will be loaded to memory
        """
        return self.n_frames * self.n_atoms * 3 * 8 / (1024 ** 3)

    @property
    def xyz(self):
        '''return 3D array of coordinates'''
        size_in_GB = self._estimated_GB
        # check if larger than size_limit_in_GB
        if size_in_GB > self._size_limit_in_GB and not self._force_load:
            raise MemoryError(
                "you are loading %s GB, larger than size_limit %s GB. "
                "Please increase self._size_limit_in_GB or set self._force_load=True"
                % (size_in_GB, self._size_limit_in_GB))
        return super(TrajectoryIterator, self).xyz

    def _iterator_slice(self, start=0, stop=None, stride=None):
        """iterator slice"""
        from itertools import islice
        if stop is None:
            stop = self.n_frames
        return islice(self, start, stop, stride)

    def _make_independent_iterators(self, n_iters):
        from itertools import tee
        return tee(self, n_iters)

    def iterframe(self,
                  start=0,
                  stop=None,
                  stride=1,
                  mask=None,
                  autoimage=False,
                  rmsfit=None,
                  copy=False,
                  frame_indices=None):
        ''''''

        if mask is None:
            _top = self.top
        else:
            _top = self.top._get_new_from_mask(mask)

        if rmsfit is not None:
            if isinstance(rmsfit, tuple):
                assert len(rmsfit) <= 2, (
                    "rmsfit must be a tuple of one (frame,) "
                    "or two elements (frame, mask)")
                if len(rmsfit) == 1:
                    rmsfit = (rmsfit[0], '*')
            elif isinstance(rmsfit, int) or isinstance(rmsfit, Frame):
                rmsfit = (rmsfit, '*')
            else:
                raise ValueError("rmsfit must be a tuple or an integer")

            if is_int(rmsfit[0]):
                index = rmsfit[0]
                rmsfit = ([self[index], rmsfit[1]])

        # check how many frames will be calculated
        if frame_indices is None:
            # only check if does not have frame_indices
            if stop is None or stop >= self.n_frames:
                stop = self.n_frames
            elif stop < 0:
                stop = get_positive_idx(stop, self.n_frames)
            n_frames = len(range(start, stop, stride))
            frame_iter_super = super(
                TrajectoryIterator, self).iterframe(start, stop, stride)
        else:
            stop = None
            start = None
            stride = None
            try:
                 n_frames = len(frame_indices)
            except TypeError:
                n_frames = None
            frame_iter_super = super(TrajectoryIterator,
                                     self)._iterframe_indices(frame_indices)

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
                         copy=copy,
                         frame_indices=frame_indices)

    def iterchunk(self,
                  chunksize=2,
                  start=0,
                  stop=-1,
                  autoimage=False,
                  rmsfit=None):
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
            for chunk in trajiter.iterchunk(100, autoimage=True, rmsfit=(ref0, '@CA'))

        Notes
        -----
        if using 'autoimage` with reference frame for rms-fit, make sure to `autoimage`
        ref first
        """
        if rmsfit is not None:
            ref, mask_for_rmsfit = rmsfit
            need_align = True
        else:
            need_align = False
            ref, mask_for_rmsfit = None, None

        for chunk in super(TrajectoryIterator, self).iterchunk(
            chunksize, start, stop):
            # always perform autoimage before doing fitting
            # chunk is `Trajectory` object, having very fast `autoimage` and
            # `rmsfit` methods
            if autoimage:
                chunk.autoimage()
            if need_align:
                chunk.superpose(ref=ref, mask=mask_for_rmsfit)
            # free memory
            # if not, memory will be quicly accumulated.
            if self._chunk:
                self._chunk.__del__()

            self._chunk = chunk
            yield self._chunk

    @property
    def filename(self):
        '''return 1st filename in filelist. For testing only
        '''
        return self.filelist[0]

    @property
    def shape(self):
        '''(n_frames, n_atoms, 3)'''
        return (self.n_frames, self.n_atoms, 3)

    def _split_iterators(self,
                        n_chunks=1,
                        start=0,
                        stop=-1,
                        stride=1,
                        mask=None,
                        autoimage=False,
                        rmsfit=None,
                        rank=0):
        """simple splitting `self` to n_chunks FrameIter objects

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> list(traj.split_iterators(n_chunks=4, mask='@CA'))
        """

        assert 0 <= start <= self.n_frames, "0 <= start <= self.n_frames"

        if stop <= 0 or stop > self.n_frames:
            stop = self.n_frames

        if n_chunks == 1:
            return self(start=start,
                       stop=stop,
                       stride=stride,
                       mask=mask,
                       autoimage=autoimage,
                       rmsfit=rmsfit)
        else:
            if rank >= 0:
                _start, _stop = split_range(n_chunks=n_chunks,
                                                   start=start,
                                                   stop=stop)[rank]
                return self.iterframe(start=_start,
                                      stop=_stop,
                                      stride=stride,
                                      mask=mask,
                                      autoimage=autoimage,
                                      rmsfit=rmsfit)
            else:
                list_of_iterators = []
                for (_start, _stop) in split_range(n_chunks=n_chunks,
                                            start=start,
                                            stop=stop):
                    list_of_iterators.append(self.iterframe(start=_start,
                                      stop=_stop,
                                      stride=stride,
                                      mask=mask,
                                      autoimage=autoimage,
                                      rmsfit=rmsfit))
                return list_of_iterators

    @property
    def temperatures(self):
        return np.array([frame.temperature for frame in self])

    def _iselect(self, frame_indices):
        return self._iterframe_indices(frame_indices)
