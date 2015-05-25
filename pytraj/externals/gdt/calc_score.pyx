from __future__ import absolute_import 
from cpython.array cimport array as pyarray
from .gdt cimport gdt
from ... _get_common_objects import _get_top, _get_data_from_dtype
from ... _shared_methods import _frame_iter_master
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
    """
    cdef pyarray score_arr = pyarray('d', [])

    if score == 'gdtscore':
        _score = 1
    elif score == 'tmscore':
        _score = 2
    elif score == 'maxsub':
        _score = 3

    if ref is None:
        _ref = traj[0]
    else:
        _ref = ref
    
    _top = _get_top(traj, top)
    atm = top(mask)
    ref_coords_1d = Frame(ref, atm).coords
    n_atoms_new = <int> len(ref_coords_1d) / 3 

    for frame in _frame_iter_master(traj):
        frame_coords_1d = Frame(frame, atm).coords
        value = gdt(ref_coords_1d, frame_coords_1d, 1, n_atoms_new, _score)[0]/1000.
        score_arr.append(value)

    if dtype == 'pyarray':
        return score_arr
    else:
        from pytraj.datasets import DataSet_double
        dset = DataSet_double()
        dset.resize(score_arr.__len__())
        dset.values[:] = score_arr
        return _get_data_from_dtype(dset, dtype)
