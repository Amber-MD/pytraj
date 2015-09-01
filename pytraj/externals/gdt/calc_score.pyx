from __future__ import absolute_import 
from cpython.array cimport array as pyarray
from .gdt cimport _gdt, sshort
from ... _get_common_objects import _get_top, _get_data_from_dtype
from ... _shared_methods import iterframe_master
from ... Frame cimport Frame

def calc_score(traj, top=None, ref=None, mask="*", 
               dtype='pyarray',
               score="gdtscore", *args, **kwd):
    """return `gdtscore` or `tmscore` or `maxsub`

    Parameters
    ---------
    traj : Trajectory-like | list-like | Frame-like
    ref : Frame, reference, default=None (use first frame)
        ref must have the same atom number as traj.top (or top)
    top : Topology object
    mask : str, atom mask
    score : str
        `gdtscore` or `tmscore` or `maxsub`

    Examples
    -------
    >>> from pytraj.common_actions import calc_protein_score
    >>> calc_protein_score(traj, ref=traj[0])
    """
    cdef pyarray score_arr = pyarray('d', [])
    cdef Frame _frame, _ref
    cdef sshort result
    cdef int int_score


    if score == 'gdtscore':
        int_score = 1
    elif score == 'tmscore':
        int_score = 2
    elif score == 'maxsub':
        int_score = 3

    if ref is None:
        try:
            _ref = traj[0]
        except IndexError:
            raise IndexError("must be a Trajectory")
    else:
        _ref = ref

    _top = _get_top(traj, top)
    atm = top(mask)
    _ref = Frame(ref, atm)
    n_atoms_new = _ref.n_atoms

    for frame in iterframe_master(traj):
        _frame = Frame(frame, atm)
        result = _gdt(_ref.thisptr.xAddress(), _frame.thisptr.xAddress(), 
                     1, n_atoms_new, int_score)[0]
        score_arr.append(<double> result/1000.)

    if dtype == 'pyarray':
        return score_arr
    else:
        from pytraj.datasets import DatasetDouble
        dset = DatasetDouble()
        dset.resize(score_arr.__len__())
        dset.values[:] = score_arr
        return _get_data_from_dtype(dset, dtype)
