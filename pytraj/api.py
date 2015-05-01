from __future__ import absolute_import
import numpy as np
from .Frame import Frame
from .Topology import Topology
from ._action_in_traj import ActionInTraj
from ._shared_methods import _frame_iter, _frame_iter_master
from .trajs.Trajin_Single import Trajin_Single
from  .utils.check_and_assert import is_int, is_frame_iter
from .externals.six import string_types
from .externals.six.moves import range
from .AtomMask import AtomMask
from ._get_common_objects  import _get_top # need to move this method to more correct module

# TODO : more checking.
class Trajectory(ActionInTraj):
    def __init__(self, filename_or_traj=None, top=None):
        self.top = _get_top(filename_or_traj, top)
        self.xyz = None

        if filename_or_traj is None or filename_or_traj == "":
            self.xyz = None
        elif hasattr(filename_or_traj, 'xyz'):
            self.xyz = filename_or_traj.xyz
        elif isinstance(filename_or_traj, string_types):
            self.load(filename_or_traj)
        else:
            raise NotImplementedError("need to have filename or 3D array or Trajectory-like object")

    def __str__(self):
        clsname = self.__class__.__name__
        txt = "%s with %s frames, %s atoms" % (clsname, self.n_frames, 
                                               self.n_atoms) 
        return txt

    def __repr__(self):
        return self.__str__()

    def copy(self):
        traj = Trajectory()
        traj.top = self.top.copy()
        traj.xyz = self.xyz.copy()
        return traj

    @property
    def shape(self):
        try:
            return self.xyz.shape
        except:
            return (None, None, 3)

    @property
    def n_atoms(self):
        try:
            n_atoms = self.xyz.shape[1]
        except:
            try:
                n_atoms = self.top.n_atoms
            except:
                n_atoms = 0
        return n_atoms

    @property
    def ndim(self):
        return 3

    @property
    def n_frames(self):
        try:
            n_frames = self.xyz.shape[0]
        except:
            n_frames = 0
        return n_frames

    def __iter__(self):
        """return a copy of Frame object"""
        for i in range(self.xyz.shape[0]):
            frame = Frame()
            frame.append_xyz(self.xyz[i])
            yield frame

    def __getitem__(self, idx):
        """return a copy of Frame object"""
        if is_int(idx):
            arr0 = self.xyz[idx]
            frame = Frame()
            frame.append_xyz(arr0)
            return frame
        else:
            traj = self.__class__()
            atm = None
            arr0 = None

            if isinstance(idx, string_types):
                atm = self.top(idx)
            elif isinstance(idx, AtomMask):
                atm = idx
            if isinstance(atm, AtomMask):
                traj.top = self.top._modify_state_by_mask(atm)
                arr0 = self.xyz[:, atm.indices]
            else:
                traj.top = self.top
                arr0 = self.xyz[idx]
            traj.append(arr0)
            return traj

    def __setitem__(self, idx, other_frame):
        """set Frame"""
        if is_int(idx):
            self.xyz[idx] = other_frame.xyz
        else:
            raise NotImplementedError("idx must be an integer")

    def append(self, other):
        """other: xyz, Frame, FrameArray, ...

        Notes
        ----
        Can not append TrajReadOnly object since we use Trajectory in TrajReadOnly class
        """
        if isinstance(other, Frame):
            arr0 = other.xyz.reshape((1, other.n_atoms, 3))
            if self.xyz is None:
                self.xyz = arr0.copy()
            else:
                self.xyz = np.append(self.xyz, arr0, axis=0)
        elif isinstance(other, np.ndarray) and other.ndim == 3:
            if self.xyz is None:
                self.xyz = other
            else:
                self.xyz = np.append(self.xyz, other, axis=0)
        elif hasattr(other, 'n_frames') and hasattr(other, 'xyz'):
            # assume Trajectory-like object
            if self.xyz is None:
                self.xyz = other.xyz
            else:
                self.xyz = np.append(self.xyz, other.xyz, axis=0)
        elif is_frame_iter(other):
            for frame in other:
                self.append(frame)
        else:
            # try to iterate to get Frame
            for frame in _frame_iter_master(other):
                self.append(frame)

    def join(self, other):
        self.append(other)

    def frame_iter(self, start=0, stop=-1, stride=1, mask=None):
        return _frame_iter(self, start, stop, stride, mask)

    def load(self, filename=''):
        ts = Trajin_Single(filename, self.top)
        self.append(ts.xyz)
