from __future__ import absolute_import
from .utils import _import_numpy
from .core import Box
from .Frame import Frame
from ._action_in_traj import ActionTrajectory
from ._shared_methods import _frame_iter, _frame_iter_master
from .trajs.TrajectoryCpptraj import TrajectoryCpptraj
from .utils.check_and_assert import is_int, is_frame_iter
from .externals.six import string_types
from .externals.six.moves import range
from .AtomMask import AtomMask
# need to move this method to more correct module
from ._get_common_objects import _get_top

_, np = _import_numpy()

__all__ = ['Trajectory']


class Trajectory(ActionTrajectory):

    def __init__(self, filename_or_traj=None, top=None):
        self.top = _get_top(filename_or_traj, top)
        self.xyz = None
        self._boxes = None
        self._life_holder = None

        if filename_or_traj is None or filename_or_traj == "":
            self.xyz = None
        elif hasattr(filename_or_traj, 'xyz'):
            # make sure to use `float64`
            self.xyz = filename_or_traj.xyz.astype(np.float64)
        elif isinstance(filename_or_traj, string_types):
            self.load(filename_or_traj)
        else:
            raise NotImplementedError(
                "need to have filename or 3D array or Trajectory-like object")
        if hasattr(filename_or_traj, 'box_to_ndarray'):
            self._boxes = filename_or_traj.box_to_ndarray()

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
            frame = Frame(self.n_atoms)
            frame[:] = self.xyz[i]
            if self._boxes is not None:
                frame.box = Box(self._boxes[i])
            self._life_holder = frame
            yield self._life_holder

    def __getitem__(self, idx):
        """return a copy of Frame object"""
        if is_int(idx):
            arr0 = self.xyz[idx]
            frame = Frame(self.n_atoms)
            frame[:] = arr0
            if self._boxes is not None:
                frame.box = Box(self._boxes[idx])
            self._life_holder = frame
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
            traj.update_box(self._boxes[idx])
            self._life_holder = traj
        return self._life_holder

    def __setitem__(self, idx, other_frame):
        """set Frame"""
        if is_int(idx):
            self.xyz[idx] = other_frame.xyz
        else:
            raise NotImplementedError("idx must be an integer")

    def __iadd__(self, other):
        self.xyz = np.vstack((self.xyz, other.xyz))
        return self

    def __add__(self, other):
        return self.__iadd__(other)

    def append(self, other):
        """other: xyz, Frame, Trajectory, ...

        Notes
        ----
        Can not append TrajectoryIterator object since we use Trajectory in TrajectoryIterator class
        """
        if isinstance(other, Frame):
            arr0 = other.xyz.reshape((1, other.n_atoms, 3))
            barr = other.box.to_ndarray().reshape((1, 6))
            if self.xyz is None:
                self.xyz = arr0.copy()
                self._boxes = barr
            else:
                self.xyz = np.vstack((self.xyz, arr0))
                self._boxes = np.vstack((self._boxes, barr))
        elif isinstance(other, np.ndarray) and other.ndim == 3:
            if self.xyz is None:
                self.xyz = other
                self._boxes = np.empty((other.shape[0], 6))
            else:
                self.xyz = np.vstack((self.xyz, other))

                fake_box_arr = np.empty((other.shape[0], 6))
                self._boxes = np.vstack((self._boxes, fake_box_arr))
        elif hasattr(other, 'n_frames') and hasattr(other, 'xyz'):
            # assume Trajectory-like object
            if self.xyz is None:
                self.xyz = other.xyz[:]
                self._boxes = other.box_to_ndarray()
            else:
                self.xyz = np.vstack((self.xyz, other.xyz))
                self._boxes = np.vstack((self._boxes, other.box_to_ndarray()))
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
        ts = TrajectoryCpptraj()
        ts.top = self.top
        ts.load(filename, top=ts.top)
        self.append(ts.xyz[:])
        self._boxes = ts.box_to_ndarray()

    def has_box(self):
        try:
            return self.top.has_box()
        except:
            return False

    def autoimage(self):
        import pytraj.common_actions as pyca
        if not self.has_box():
            print("there is no box, skip")
        else:
            for idx, frame in enumerate(self):
                pyca.autoimage(frame, top=self.top)
                self.xyz[idx] = frame.xyz[:]

    def box_to_ndarray(self):
        return self._boxes

    def update_box(self, box_arr):
        self._boxes = box_arr
