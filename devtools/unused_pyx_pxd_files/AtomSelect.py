"""mimic Atom Select in VMD"""
from __future__ import absolute_import
from pytraj.externals.six import string_types
from array import array as pyarray
#from .Frame import Frame
from .utils.check_and_assert import _import_numpy
from .AtomMask import AtomMask

# TODO : using 'vmd_style' flag to turn on/off mask style

has_numpy, np = _import_numpy()
if not has_numpy:
    raise ImportError("must has numpy installed")

# make sure to inherit from `object` class so we can use property


class AtomSelect(object):

    def __init__(self, top=None, traj=None, frameidx=0, frame=None):
        self.top = top
        self.traj = traj
        self._frameidx = frameidx
        if self.traj is not None:
            self._selected_frame = self.traj[self._frameidx]

    def set_top(self, top):
        # we use "set" for those methods to emphasize
        self.top = top

    @property
    def frameidx(self):
        return self._frameidx

    @frameidx.setter
    def frameidx(self, value):
        self._frameidx = value
        # set self._selected_frame too
        self._selected_frame = self.traj[value]

    def set_traj(self, traj):
        self.traj = traj

    @property
    def selected_frame(self):
        return self._selected_frame

    @selected_frame.setter
    def selected_frame(self, value):
        # FIXME: update self.frameidx too
        self._selected_frame = value

    def __call__(self, *args, **kwd):
        return self.select(*args, **kwd)

    def select(self, mask, frameidx=None):
        """return 2D numpy array"""
        if frameidx is not None:
            self.frameidx = frameidx
        arr0 = np.asarray(self._selected_frame.buffer2d)

        if isinstance(mask, AtomMask):
            return arr0[mask.selected_indices()]
        elif isinstance(mask, string_types):
            atm = AtomMask(mask)
            self.top.set_integer_mask(atm)
            return arr0[atm.selected_indices()]
        elif isinstance(mask, (pyarray, np.ndarray, list, tuple)):
            return arr0[mask]
        else:
            raise NotImplementedError(
                "must be string, AtomMask or arrayl-like")

    def get_indices(self, mask):
        """return a list of atom indinces"""

        if isinstance(mask, AtomMask):
            return mask.selected_indices()
        elif isinstance(mask, string_types):
            atm = AtomMask(mask)
            self.top.set_integer_mask(atm)
            return atm.selected_indices()
