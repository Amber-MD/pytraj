"""out-of-core TrajectoryIterator
"""
from __future__ import absolute_import
import os
import re
from glob import glob
import numpy as np
from .trajs.cpptraj_trajectory import TrajectoryCpptraj
from .externals.six import string_types
from .externals.six.moves import range
from .topology import Topology
from .frame import Frame
from .utils import is_int
from .cyutils import get_positive_idx
from .frameiter import FrameIterator
from ._get_common_objects import _load_Topology
from .utils import split_range
from .utils.convert import array_to_cpptraj_atommask

__all__ = ['TrajectoryIterator', ]


# tryint, alphanum_key, sort_filename_by_number are adapted from
# http://nedbatchelder.com/blog/200712/human_sorting.html
def tryint(s):
    try:
        return int(s)
    except:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [tryint(c) for c in re.split('([0-9]+)', s)]


def sort_filename_by_number(filelist):
    """ Sort the given list in the way that humans expect.
    """
    return sorted(filelist, key=alphanum_key)


def _make_frame_slices(n_files, original_frame_slice):
    '''
    >>> _make_frame_slices(2, (0, -1))
    [(0, -1), (0, -1)]

    >>> _make_frame_slices(3, [(0, -1), (0, -2)],)
    [(0, -1), (0, -2), (0, -1, 1)]

    >> # raise
    >>> _make_frame_slices(3, None)
    Traceback (most recent call last):
        ...
    ValueError: must be a tuple of integer values or a list of tuple of integer values

    >>> _make_frame_slices(2, [(0, -1), (0, -2), (0, -1, 1)])
    Traceback (most recent call last):
        ...
    ValueError: len of frame_slice tuple-list must be smaller or equal number of files
    '''
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
    '''out-of-core trajectory holder.

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.testing import get_fn
    >>> traj_name, top_name = get_fn('tz2')
    >>> traj = pt.TrajectoryIterator(traj_name, top_name)

    >>> # user should always use :method:`pytraj.iterload` to load TrajectoryIterator
    >>> traj = pt.iterload(['remd.x.000', 'remd.x.001'], 'input.parm7') # doctest: +SKIP

    Notes
    -----
    It's a bit tricky to pickle this class. As default, new TrajectoryIterator will
    use original trj filename and top filename. If set _pickle_topology to True, its
    Topology will be pickled (slow but more accurate if you change the topology in the
    fly)
    '''

    def __init__(self, filename=None, top=None, *args, **kwd):
        self._force_load = False
        # use self._chunk to store `chunk` in iterchunk
        # to deallocate memory
        self._chunk = None
        # only allow to load <= 1 GB
        self._size_limit_in_GB = 1
        self._pickle_topology = False
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
            self._load(filename, self.top, *args, **kwd)
        if not top and (args or kwd):
            raise ValueError('require a Topology')

        self.__dict__.update({
            '_top_filename': self.top.filename,
            'filelist': self.filelist,
            'frame_slice_list': self.frame_slice_list,
        })

    def __setstate__(self, state):
        self.__dict__ = state
        if self._pickle_topology:
            self.top = state['top']
        else:
            # faster
            self.top = _load_Topology(state['_top_filename'])
        self._load(state['filelist'], frame_slice=state['frame_slice_list'])

    def __getstate__(self):
        '''

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_ala3()
        >>> # pickle by reloading Topology from filename
        >>> pt.to_pickle(traj, 'output/test.pk')

        >>> # pickle by rebuilding Topology
        >>> traj._pickle_topology = True
        >>> pt.to_pickle(traj, 'output/test.pk')
        '''
        if 'top' not in self.__dict__.keys() and self._pickle_topology:
            # slow
            # Topology is pickable
            if self._pickle_topology:
                self.__dict__.update({'top': self.top})

        return self.__dict__

    def __iter__(self):
        '''do not make a frame copy here
        '''
        for frame in super(TrajectoryIterator, self).__iter__():
            yield frame

    def copy(self):
        '''return a deep copy. Use this method with care since the copied traj just reuse
        the filenames
        '''
        other = self.__class__()
        other.top = self.top.copy()

        for fname, frame_slice in zip(self.filelist, self.frame_slice_list):
            other._load(fname, frame_slice=frame_slice)
        return other

    def _load(self, filename=None, top=None, frame_slice=(0, -1, 1)):
        """load trajectory/trajectories from filename/filenames
        with a single frame_slice or a list of frame_slice
        """
        if not top:
            _top = self.top
        else:
            _top = top

        if isinstance(filename, string_types) and os.path.exists(filename):
            super(TrajectoryIterator, self)._load(filename, _top, frame_slice)
            self.frame_slice_list.append(frame_slice)
        elif isinstance(filename,
                        string_types) and not os.path.exists(filename):
            flist = sort_filename_by_number(glob(filename))
            if not flist:
                raise ValueError(
                    "must provie a filename or list of filenames or file pattern")
            self._load(flist, top=top, frame_slice=frame_slice)
        elif isinstance(filename, (list, tuple)):
            filename_list = filename
            full_frame_slice = _make_frame_slices(
                len(filename_list), frame_slice)

            for fname, fslice in zip(filename_list, full_frame_slice):
                self.frame_slice_list.append(frame_slice)
                super(TrajectoryIterator, self)._load(fname,
                                                      _top,
                                                      frame_slice=fslice)
        else:
            raise ValueError("filename must a string or a list of string")

    @property
    def topology(self):
        """traditional name for Topology file

        Examples
        --------
        >>> import pytraj as pt
        >>> from pytraj.testing import get_fn
        >>> fname, tname = get_fn('ala3')
        >>> traj = pt.iterload(fname, tname)
        >>> traj.topology
        <Topology: 34 atoms, 3 residues, 1 mols, non-PBC>
        >>> new_traj = pt.TrajectoryIterator()
        >>> new_traj.topology = traj.topology
        >>> new_traj._load(traj.filename)
        """
        return self.top

    @topology.setter
    def topology(self, newtop):
        self.top = newtop

    @property
    def _estimated_GB(self):
        """esimated GB of data will be loaded to memory
        """
        return self.n_frames * self.n_atoms * 3 * 8 / (1024**3)

    @property
    def xyz(self):
        '''return 3D array of coordinates'''
        size_in_GB = self._estimated_GB
        # check if larger than size_limit_in_GB
        if size_in_GB > self._size_limit_in_GB and not self._force_load:
            raise MemoryError(
                "you are loading {0} GB, larger than size_limit {1} GB. "
                "Please increase traj._size_limit_in_GB or set traj._force_load to True".format(size_in_GB, self._size_limit_in_GB))
        return super(TrajectoryIterator, self).xyz

    def iterframe(self,
                  start=0,
                  stop=None,
                  step=1,
                  mask=None,
                  autoimage=False,
                  rmsfit=None,
                  copy=False,
                  frame_indices=None):
        '''iterate trajectory with given frame_indices or given (start, stop, step)

        Parameters
        ----------
        start : int, default 0
        stop : {None, int}, default None
            if None, iterate to final frame
        step : int, default 1
        mask : {None, str}, default None
            if None, use all atoms. If not None, use given mask
        autoimage : bool, default False
            if True, perform autoimage for each frame
        rmsfit : {None, int, tuple}, default None
            if not None, perform superpose each Frame to to reference.
        frame_indices : {None, array-like}
            if not None, iterate trajectory for given indices. If frame_indices is given,
            (start, stop, step) will be ignored.

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> for frame in traj.iterframe(0, 8, 2): pass
        >>> for frame in traj.iterframe(0, 8, 2, autoimage=True): pass
        >>> # use negative index
        >>> traj.n_frames
        10
        >>> fi = traj.iterframe(0, -1, 2, autoimage=True)
        >>> fi.n_frames
        5
        >>> # mask is atom indices
        >>> fi = traj.iterframe(0, -1, 2, mask=range(100), autoimage=True)
        >>> fi.n_atoms
        100
        '''

        if mask is None:
            _top = self.top
        else:
            if isinstance(mask, string_types):
                mask = mask
                _top = self.top._get_new_from_mask(mask)
            else:
                mask = array_to_cpptraj_atommask(mask)
                _top = self.top._get_new_from_mask(mask)

        if rmsfit is not None:
            if isinstance(rmsfit, tuple):
                assert len(rmsfit) == 2, (
                    "rmsfit must be a tuple of one (frame,) "
                    "or two elements (frame, mask)")
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
            n_frames = len(range(start, stop, step))
            frame_iter_super = super(TrajectoryIterator,
                                     self).iterframe(start, stop, step)
        else:
            stop = None
            start = None
            step = None
            try:
                n_frames = len(frame_indices)
            except TypeError:
                # itertools.chain
                n_frames = None
            frame_iter_super = super(TrajectoryIterator,
                                     self)._iterframe_indices(frame_indices)

        return FrameIterator(frame_iter_super,
                             original_top=self.top,
                             new_top=_top,
                             start=start,
                             stop=stop,
                             step=step,
                             mask=mask,
                             autoimage=autoimage,
                             rmsfit=rmsfit,
                             n_frames=n_frames,
                             copy=copy,
                             frame_indices=frame_indices)

    def iterchunk(self,
                  chunksize=2,
                  start=0,
                  stop=-1,
                  autoimage=False,
                  rmsfit=None):
        """iterate trajectory by chunk

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
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> ref = traj[3]
        >>> for chunk in traj.iterchunk(3, autoimage=True, rmsfit=(ref, '@CA')): pass

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

        for chunk in super(TrajectoryIterator, self).iterchunk(chunksize,
                                                               start, stop):
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
        '''(n_frames, n_atoms, 3)

        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2_ortho()
        >>> traj.shape
        (10, 5293, 3)
        '''
        return (self.n_frames, self.n_atoms, 3)

    def _split_iterators(self,
                         n_chunks=1,
                         start=0,
                         stop=-1,
                         step=1,
                         mask=None,
                         autoimage=False,
                         rmsfit=None,
                         rank=0):
        """simple splitting `self` to n_chunks FrameIterator objects

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> list(traj._split_iterators(n_chunks=4, mask='@CA'))
        [<Frame with 12 atoms>, <Frame with 12 atoms>]
        >>> isinstance(traj._split_iterators(n_chunks=4, mask='@CA', rank=-1), list)
        True

        >>> # reset stop value to max n_framaes if this number looks 'weird'
        >>> fi = traj._split_iterators(n_chunks=4, mask='@CA', stop=-100)
        """

        assert 0 <= start <= self.n_frames, "0 <= start <= self.n_frames"

        if stop <= 0 or stop > self.n_frames:
            stop = self.n_frames

        if rank >= 0:
            _start, _stop = split_range(n_chunks=n_chunks,
                                        start=start,
                                        stop=stop)[rank]
            return self.iterframe(start=_start,
                                  stop=_stop,
                                  step=step,
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
                                                        step=step,
                                                        mask=mask,
                                                        autoimage=autoimage,
                                                        rmsfit=rmsfit))
            return list_of_iterators

    @property
    def temperatures(self):
        return np.array([frame.temperature for frame in self])
