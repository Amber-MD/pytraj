# distutils: language = c++
#
cimport cython
from pytraj.Frame cimport _Frame, Frame
from pytraj.AtomMask cimport AtomMask
from pytraj.trajs.Trajout import Trajout
from pytraj.externals.six import string_types
from pytraj.six_2 import set
from pytraj.utils import _import_numpy
from pytraj.exceptions import PytrajMemviewError, PytrajConvertError
from pytraj.utils.check_and_assert import is_frame_iter

def _savetraj(self, filename="", fmt='unknown', overwrite=False):
    if fmt == 'unknown':
        # convert to "UNKNOWN_TRAJ"
        fmt = fmt.upper() + "_TRAJ"
    else:
        fmt = fmt.upper()

    with Trajout(filename=filename, top=self.top, fmt=fmt, 
                 overwrite=overwrite, more_args=None) as trajout:
        for idx, frame in enumerate(self):
            trajout.writeframe(idx, frame, self.top)

def _get_temperature_set(self):
    return set(self.temperatures) 

def _xyz(self):
    has_np, np = _import_numpy()
    if has_np:
        return self[:, :, :]
    else:
        raise NotImplementedError("require numpy")

def _tolist(self):
    """return flatten list for traj-like object"""
    from itertools import chain
    return [frame.tolist() for frame in self]

def my_str_method(self):
    name = self.__class__.__name__
    tmps = """%s with %s frames, %s atoms/frame
           """ % (
            name, self.size, self.top.n_atoms,
            )
    return tmps

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.infer_types(True)
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

    obj : could be FrameArray, TrajReadOnly, TrajinList, frame_iter object
          could be a list of trajs, ...
    """
    cdef Frame frame
    cdef object traj_obj

    is_frame_iter_but_not_master = (is_frame_iter(obj) and not obj.__name__ is '_frame_iter_master')
    if hasattr(obj, 'n_frames') or is_frame_iter_but_not_master:
        # traj-like or frame_iter or _frame_iter
        for frame in obj:
            yield frame
    else:
        try:
            for traj_obj in obj:
                for frame in traj_obj:
                    yield frame
        except:
            raise PytrajConvertError("can not convert to Frame")

