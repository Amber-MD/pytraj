from __future__ import absolute_import

import numpy as np
from .core.Box import Box
from .Frame import Frame
from .utils.check_and_assert import is_int, is_frame_iter
from .externals.six import string_types
from .externals.six.moves import range
from .core.cpptraj_core import AtomMask

# use absolute import here
from pytraj._get_common_objects import _get_top

from .Topology import Topology
from ._shared_methods import _savetraj, iterframe_master, my_str_method
from ._fast_iterframe import _fast_iterptr, _fast_iterptr_withbox
from .frameiter import FrameIter

__all__ = ['Trajectory']


class Trajectory(object):
    def __init__(self, filename=None, top=None, xyz=None, indices=None):
        """
        Notes
        -----
        This class is not very well-tested (comparead to pytraj.Trajectory and 
        pytraj.TrajectoryIterator classes).

        Use it with your own risk. It's just a simple numpy wrapping (like mdtraj).
        I myself do not use this much since it requires extra copying to Frame class.
        However, if you like numpy and fancy indexing, you should definitely use this.
        This is good for any cpptraj calculation without changing coordinates.

        For `pytraj.Trajectory` and `pytraj.TrajectoryIterator`, any time calling `xyz` attribute,
        a copy of xyz coordinates will be returned but this class return a memmoryview.
        
        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.Trajectory(xyz, top)
        >>> traj = pt.Trajectory(traj0, top)
        >>> traj = pt.Trajectory("traj.x", "t.prmtop"))

        >>> traj['@CA'].xyz[:, :, 0]

        """
        self._top = _get_top(filename, top)

        if self._top is None:
            self._top = Topology()

        self._xyz = None
        self._boxes = None

        # use those to keep lifetime of Frame
        self._life_holder = None
        self._frame_holder = None

        if filename is None or filename == "":
            if xyz is not None:
                if self.top.is_empty():
                    raise ValueError("must have a non-empty Topology")
                else:
                    assert self.top.n_atoms == xyz.shape[
                        1
                    ], "must have the same n_atoms"
                self._xyz = np.asarray(xyz)
            else:
                self._xyz = None
        elif hasattr(filename, 'xyz'):
            # make sure to use `float64`
            self._xyz = filename.xyz.astype(np.float64)
        elif isinstance(filename, (string_types, list, tuple)):
            if isinstance(filename, string_types):
                self.load(filename)
            else:
                for fname in filename:
                    self.load(fname)

        elif is_frame_iter(filename):
            for frame in filename:
                self.append(frame.xyz[:])
        else:
            self._xyz = np.asarray(filename, dtype='f8')

        if hasattr(self._xyz, 'shape'):
            assert self.top.n_atoms == self._xyz.shape[
                1
            ], "must have the same n_atoms"

        if hasattr(filename, 'unitcells'):
            self._boxes = filename.unitcells

    @property
    def top(self):
        return self._top

    @top.setter
    def top(self, value):
        self._top = value.copy()

    def reverse(self):
        self._xyz = self._xyz[::-1]
        if self._boxes is not None:
            self._boxes = self._boxes[::-1]

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, values):
        if self.shape[1]:
            if self.n_atoms != values.shape[1]:
                raise ValueError("must have the same number of atoms")
        if not values.flags['C_CONTIGUOUS']:
            raise TypeError('must be C_CONTIGUOUS')
        self._xyz = values

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return self.__str__()

    def copy(self):
        traj = self.__class__()
        traj.top = self.top.copy()
        traj.xyz = self._xyz.copy()
        return traj

    @property
    def shape(self):
        '''(n_frames, n_atoms, 3)
        '''
        try:
            return self._xyz.shape
        except:
            return (None, None, 3)

    @property
    def n_atoms(self):
        '''n_atoms
        '''
        return self.top.n_atoms

    @property
    def n_frames(self):
        '''n_frames
        '''
        try:
            n_frames = self._xyz.shape[0]
        except (AttributeError, IndexError):
            n_frames = 0
        return n_frames

    @property
    def size(self):
        '''alias of Trajectory.n_frames
        '''
        return self.n_frames

    def __iter__(self):
        """return a Frame view of coordinates

        Notes
        -----
        update frame view will update Trajectory.xyz too
        if want to use listcomp, need to make copy for every frame
        >>> [frame.copy() for frame in traj]

        Examples
        --------
        >>> for frame in traj: print(frame)
        """
        indices = range(self.n_frames)
        return self._iterframe_indices(indices)

    def _iterframe_indices(self, indices):
        """return a Frame view of coordinates

        Notes
        -----
        update frame view will update Trajectory.xyz too

        Examples
        --------
        >>> for frame in traj: print(frame)
        """

        if self._boxes is None:
            return _fast_iterptr(self.xyz, self.n_atoms, indices)
        else:
            return _fast_iterptr_withbox(self.xyz, self._boxes, self.n_atoms,
                                         indices)

    def __getitem__(self, idx):
        """return a view or copy of coordinates (follow numpy's rule)

        Examples
        --------
        >>> # create mutable trajectory from TrajectoryIterator
        >>> t0 = traj[:]
        >>> print(t0)

        >>> # get a Frame view
        >>> t0[2]

        >>> # get a Trajetory view
        >>> t0[0:8:2]

        >>> # get a copy of Trajetory
        >>> t0[[0, 4, 6]]

        >>> # get a copy, keep only CA atoms
        >>> t0['@CA']

        >>> # get a copy, keep only CA atoms for 3 frames
        >>> t0[:3, '@CA']

        >>> # get a new stripped Frame
        >>> t0[0, '@CA']
        """
        if is_int(idx):
            # traj[0]
            # return a single Frame as a view
            arr0 = self._xyz[idx]
            frame = Frame(self.n_atoms, arr0, _as_ptr=True)
            if self._boxes is not None:
                frame.box = Box(self._boxes[idx])
            self._life_holder = frame
        else:
            # return a new Trajectory
            traj = self.__class__()
            atm = None
            arr0 = None

            if isinstance(idx, (string_types, AtomMask)):
                # return a copy
                # traj['@CA']
                if isinstance(idx, string_types):
                    atm = self.top(idx)
                elif isinstance(idx, AtomMask):
                    atm = idx
                else:
                    raise IndexError("only support integer or string indexing")
                if isinstance(atm, AtomMask):
                    traj.top = self.top._modify_state_by_mask(atm)
                    arr0 = self._xyz[:, atm.indices]
                else:
                    traj.top = self.top
                    arr0 = self._xyz[idx]
                # make copy to create contigous memory block
                traj._xyz = arr0.copy()

                if self._boxes is not None:
                    # always make a copy in this case
                    traj._boxes = self._boxes.copy()
            elif not isinstance(idx, tuple):
                # might return a view or a copy
                # based on numpy array rule
                # traj.xyz[idx]
                traj.top = self.top
                traj._xyz = self._xyz[idx]
                if self._boxes is not None:
                    traj._boxes = self._boxes[idx]
            else:
                # is a tuple
                if len(idx) == 1:
                    traj = self[idx[0]]
                elif len(idx) == 2 and is_int(idx[0]) and isinstance(
                        idx[1], string_types):
                    # traj[0, '@CA']: return a stripped Frame
                    frame = self[idx[0]].copy()
                    # make AtomMask object
                    atm = self.top(idx[1])
                    # create new Frame with AtomMask
                    self._life_holder = Frame(frame, atm)
                    return self._life_holder
                else:
                    self._life_holder = self[idx[0]]
                    if isinstance(self._life_holder, Frame):
                        self._frame_holder = self._life_holder
                    traj = self._life_holder[idx[1:]]
            self._life_holder = traj
        return self._life_holder

    def __setitem__(self, idx, other):
        if self.n_frames == 0:
            raise ValueError("Your Trajectory is empty, how can I index it?")

        if other is None:
            raise ValueError("why bothering assign None?")
        if is_int(idx):
            if hasattr(other, 'xyz') or isinstance(other, Frame):
                # traj[1] = frame
                self._xyz[idx] = other.xyz
            else:
                # traj[1] = xyz
                # check shape?
                self._xyz[idx] = other
        elif idx == '*':
            # why need this?
            # traj.xyz = xyz
            # update all atoms, use fast version
            self._xyz[:] = other  # xyz
        elif isinstance(idx, AtomMask) or isinstance(idx, string_types):
            # update xyz for mask
            # traj['@CA'] = xyz
            if isinstance(idx, AtomMask):
                atm = idx
            else:
                atm = self.top(idx)
            if isinstance(other, Trajectory):
                indices = atm.indices

                for i in range(self.n_frames):
                    for j, k in enumerate(indices):
                        self.xyz[i, k] = other.xyz[i, j]
            else:
                view3d = other
                try:
                    int_view = atm.indices.astype('i4')
                except ValueError:
                    int_view = atm.indices
                # loop all frames
                for i in range(view3d.shape[0]):
                    self._xyz[:, int_view] = view3d[:]
        else:
            # really need this?
            # example: self[0, 0, 0] = 100.
            self._xyz[idx] = other

    def __iadd__(self, other):
        '''use with care. not for regular user
        '''
        if hasattr(other, 'xyz'):
            self._xyz.__iadd__(other.xyz)
        else:
            self._xyz.__iadd__(other)
        return self

    def __add__(self, other):
        '''use with care. not for regular user
        '''
        self_cp = self.copy()
        self_cp.__iadd__(other)
        return self_cp

    def append_xyz(self, xyz):
        '''append 3D numpy array
        '''
        # make sure 3D
        if xyz.ndim != 3:
            raise ValueError("ndim must be 3")

        if self.shape == (None, None, 3):
            self._xyz = xyz
        else:
            self._xyz = np.vstack((self._xyz, xyz))

    def _append_unitcells(self, box):
        if isinstance(box, tuple):
            clen, cangle = box
            data = np.hstack((clen, cangle))
            if self._boxes is None:
                self._boxes = np.asarray([data])
            else:
                self._boxes = np.vstack((self._boxes, data))

            if self._boxes.ndim == 3:
                self._boxes = self._boxes.reshape((self.n_frames, 6))
        else:
            self._boxes = np.vstack((self._boxes, box))

    def append(self, other):
        """other: xyz, Frame, Trajectory, ...

        Examples
        --------
        >>> f0 = traj0[0]
        >>> traj1.append(f0)

        Notes
        -----
        Can not append TrajectoryIterator object
        since we use Trajectory in TrajectoryIterator class
        """
        if isinstance(other, Frame):
            arr0 = other.xyz.reshape((1, other.n_atoms, 3))
            barr = other.box.to_ndarray().reshape((1, 6))
            if self._xyz is None:
                self._xyz = arr0.copy()
                self._boxes = barr
            else:
                self._xyz = np.vstack((self._xyz, arr0))
                self._boxes = np.vstack((self._boxes, barr))
        elif isinstance(other, np.ndarray) and other.ndim == 3:
            if self._xyz is None:
                self._xyz = other
                self._boxes = np.empty((other.shape[0], 6))
            else:
                self._xyz = np.vstack((self._xyz, other))

                fake_box_arr = np.empty((other.shape[0], 6))
                self._boxes = np.vstack((self._boxes, fake_box_arr))
        elif hasattr(other, 'n_frames') and hasattr(other, 'xyz'):
            # assume Trajectory-like object
            if self._xyz is None:
                self._xyz = other.xyz[:]
                self._boxes = other.unitcells
            else:
                self._xyz = np.vstack((self._xyz, other.xyz))
                self._boxes = np.vstack((self._boxes, other.unitcells))
        elif is_frame_iter(other):
            for frame in other:
                self.append(frame)
        else:
            # try to iterate to get Frame
            for frame in iterframe_master(other):
                self.append(frame)

    def join(self, other):
        if isinstance(other, Trajectory):
            self.append_xyz(other.xyz)
            if self.unitcells is not None and other.unitcells is not None:
                self._append_unitcells(other.unitcells)
        else:
            ValueError()

    def __call__(self, *args, **kwd):
        return self.iterframe(*args, **kwd)

    def _load_new_by_scipy(self, filename):
        from scipy import io
        import numpy as np

        fh = io.netcdf_file(filename, mmap=False)
        self.xyz = fh.variables['coordinates'].data
        cell_lengths = fh.variables['cell_lengths'].data
        cell_angles = fh.variables['cell_angles'].data
        self.unitcells = np.hstack((cell_lengths, cell_angles))

    def load(self, filename='', top=None, indices=None):
        if top is not None:
            if self.top.is_empty():
                self.top = top.copy()
            else:
                pass
            # don't update top if not self.top.is_empty()
        else:
            if self.top.is_empty():
                # if both top and self.top are empty, need to raise ValueError
                try:
                    tmpobj = filename
                    if hasattr(tmpobj, 'top'):
                        self.top = tmpobj.top.copy()
                    elif hasattr(tmpobj[0], 'top'):
                        self.top = tmpobj[0].top.copy()
                except:
                    raise ValueError("need to have non-empty Topology")

        # always use self.top
        if isinstance(filename, string_types):
            from pytraj.TrajectoryIterator import TrajectoryIterator
            ts = TrajectoryIterator()
            ts.top = self.top.copy()
            ts.load(filename)
            if indices is None:
                self.append_xyz(ts.xyz)
            elif isinstance(indices, slice):
                self.append_xyz(ts[indices].xyz)
            else:
                # indices is tuple, list, ...
                # we loop all traj frames and extract frame-ith in indices 
                # TODO : check negative indexing?
                # increase size of vector
                for idx in indices:
                    self.append_xyz(ts[idx].xyz)
        elif isinstance(filename, Frame):
            self.append(filename)
        elif isinstance(filename, (list, tuple)):
            # load from a list/tuple of filenames
            # or a list/tuple of numbers
            _f0 = filename[0]
            if isinstance(_f0, string_types) or hasattr(_f0, 'n_frames'):
                # need to check `string_types` since we need to load list of numbers too.
                # list of filenames
                list_of_files_or_trajs = filename
                for fh in list_of_files_or_trajs:
                    if self.warning:
                        print("Loading from list/tuple. Ignore `indices`")
                    # recursive
                    self.load(fh, self.top, indices)
            else:
                # load xyz
                try:
                    _xyz = filename
                    self.append_xyz(_xyz)
                except:
                    raise ValueError(
                        "must be a list/tuple of either filenames/Traj/numbers")
        elif hasattr(filename, 'n_frames'):
            # load from Traj-like object
            # make temp traj to remind about traj-like
            traj = filename
            if indices is None:
                for frame in traj:
                    self.append(frame)
            else:
                for idx, frame in enumerate(traj):
                    # slow method.
                    if idx in indices:
                        self.append(frame)
        elif is_frame_iter(filename):
            # load from frame_iter
            _frame_iter = filename
            for frame in _frame_iter:
                self.append(frame)
        else:
            try:
                # load from array
                _xyz = filename
                self.append_xyz(_xyz)
            except:
                raise ValueError(
                    "filename must be str, traj-like or numpy array")

        try:
            if self._xyz.shape != self.unitcells.shape:
                print("make sure to update traj.unitcells too")
        except AttributeError:
            print("make sure to update traj.unitcells too")

    def has_box(self):
        try:
            return self.top.has_box()
        except:
            return False

    def center(self, mask="", *args, **kwd):
        '''
        '''
        from pytraj.actions.CpptrajActions import Action_Center as Action

        act = Action()
        act.read_input(mask, top=self.top)
        act.process(self.top)

        for idx, frame in enumerate(self):
            frame.set_frame_mass(self.top)
            act.do_action(frame)
            self._xyz[idx] = frame.xyz[:]

    def autoimage(self):
        from pytraj.actions.CpptrajActions import Action_AutoImage

        if not self.has_box():
            raise ValueError("must have a box")
        else:
            act = Action_AutoImage()
            act.read_input("", top=self.top)
            act.process(self.top)

            for idx, frame in enumerate(self):
                act.do_action(frame)

    def rotate(self, *args, **kwd):
        import pytraj.common_actions as pyca

        for idx, frame in enumerate(self):
            pyca.rotate(frame, top=self.top, *args, **kwd)
            self.xyz[idx] = frame.xyz

    @property
    def unitcells(self):
        return self._boxes

    @unitcells.setter
    def unitcells(self, values):
        self._boxes = values

    def rmsfit(self, ref=None, mask="*"):
        """do the fitting to reference Frame by rotation and translation
        Parameters
        ----------
        ref : {Frame object, int, str}, default=None 
            Reference
        mask : str or AtomMask object, default='*' (fit all atoms)

        Examples
        --------
            traj.rmsfit(0) # fit to 1st frame
            traj.rmsfit('last', '@CA') # fit to last frame using @CA atoms
        """
        # not yet dealed with `mass` and box
        from pytraj.actions.CpptrajActions import Action_Rmsd

        if isinstance(ref, Frame):
            ref_frame = ref
        elif is_int(ref):
            i = ref
            ref_frame = self[i]
        elif isinstance(ref, string_types):
            if ref.lower() == 'first':
                i = 0
            if ref.lower() == 'last':
                i = -1
            ref_frame = self[i]
        else:
            # first
            ref_frame = self[0]

        if isinstance(mask, string_types):
            atm = self.top(mask)
        elif isinstance(mask, AtomMask):
            atm = mask
        else:
            raise ValueError("mask must be string or AtomMask object")

        for idx, frame in enumerate(self):
            _, mat, v1, v2 = frame.rmsd(ref_frame, atm, get_mvv=True)
            frame.trans_rot_trans(v1, mat, v2)

    def _allocate(self, n_frames, n_atoms):
        '''allocate (n_frames, n_atoms, 3) coordinates
        '''
        self._xyz = np.zeros((n_frames, n_atoms, 3), dtype='f8')

    def strip_atoms(self, mask):
        self.strip(mask)

    def strip(self, mask):
        '''strip atoms with given mask

        Examples
        --------
        >>> traj.strip('!@CA') # keep only CA atoms
        '''
        # AtomMask
        atm = self.top(mask)
        atm.invert_mask()
        self.top.strip_atoms(mask)

        if self._xyz is not None:
            # need to copy to make contigous memory block
            self._xyz = self._xyz[:, atm.indices].copy()

    def save(self,
             filename="",
             format='unknown',
             overwrite=True, *args, **kwd):
        _savetraj(self, filename, format, overwrite, *args, **kwd)

    def iterframe(self,
                  start=0,
                  stop=None,
                  stride=1,
                  mask=None,
                  autoimage=False,
                  frame_indices=None,
                  rmsfit=None,
                  copy=False):


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
                    rmsfit = (rmsfit, '*')
            elif isinstance(rmsfit, int):
                rmsfit = (rmsfit, '*')
            else:
                raise ValueError("rmsfit must be a tuple or an integer")

            if is_int(rmsfit[0]):
                index = rmsfit[0]
                rmsfit = ([self[index], rmsfit[1]])

        # check how many frames will be calculated
        if frame_indices is None:
            if stop is None or stop >= self.n_frames:
                stop = self.n_frames
            elif stop < 0:
                stop = get_positive_idx(stop, self.n_frames)
            else:
                stop = stop

            # make sure `range` return iterator
            indices = range(start, stop, stride)
            n_frames = len(indices)
        else:
            # frame_indices is not None
            start, stop, stride = None, None, None
            n_frames = len(frame_indices)
            indices = frame_indices

        frame_iter_super = self._iterframe_indices(indices)

        return FrameIter(frame_iter_super,
                         original_top=self.top,
                         new_top=_top,
                         start=start,
                         stop=stop,
                         stride=stride,
                         mask=mask,
                         autoimage=autoimage,
                         rmsfit=rmsfit,
                         n_frames=n_frames,
                         copy=copy)

    @property
    def _estimated_MB(self):
        """esimated MB of data will be loaded to memory
        """
        return self.n_frames * self.n_atoms * 3 * 8 / (1024 ** 2)

    @classmethod
    def from_iterable(cls, iterables, top=None, n_frames=None):
        '''
        >>> itertraj = pt.iterload('traj.nc', 'parm.top')
        >>> pt.Trajectory.from_iterable(itertraj(3, 8, 2))
        '''
        if top is None or top.is_empty():
            if hasattr(iterables, 'top'):
                top = iterables.top
            else:
                raise ValueError("must provide non-empty Topology")

        fa = Trajectory()
        fa.top = top

        if n_frames is not None:
            _n_frames = n_frames
        elif hasattr(iterables, 'n_frames'):
            _n_frames = iterables.n_frames
        else:
            try:
                _n_frames = len(iterables)
            except:
                _n_frames = None

        if _n_frames is None:
            for frame in iterables:
                # slow
                fa.append(frame)
        else:
            # faster
            fa._allocate(_n_frames, fa.top.n_atoms)
            fa._boxes = np.empty((_n_frames, 6), dtype='f8')
            for idx, frame in enumerate(iterables):
                fa._xyz[idx] = frame.xyz
                fa._boxes[idx] = frame.box.data
        return fa

    def __len__(self):
        return self.n_frames

    def __del__(self):
        self._xyz = None
        self._boxes = None

    def _apply(self, func):
        for x in self.xyz:
            x = func(x)
