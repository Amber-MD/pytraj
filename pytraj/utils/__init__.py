""""""
from __future__ import absolute_import
from .check_and_assert import assert_almost_equal, file_exist, is_generator
from .check_and_assert import assert_almost_equal as aa_eq
from .check_and_assert import eq
from .check_and_assert import _import, is_int
from .check_and_assert import has_, is_array
from .check_and_assert import ensure_not_none_or_string
from .Timer import Timer
from .context import tempfolder
from ..externals.six.moves import range
from . import convert


def duplicate_traj(orig_traj, n_times):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('ala3')
    >>> traj.n_frames
    1
    >>> duplicate_traj(traj, 3).n_frames
    3
    >>> t0 = pt.Trajectory(xyz=traj.xyz, top=traj.top)
    >>> duplicate_traj(t0, 3).n_frames
    4
    '''
    traj = orig_traj.copy()
    for _ in range(n_times - 1):
        if 'Iter' in orig_traj.__class__.__name__:
            # TrajectoryIterator
            traj._load(orig_traj.filelist)
        else:
            # Trajectory
            traj.append(traj.copy())
    return traj


def join_mask(m, res=None):
    """
    Examples
    --------
    >>> join_mask(('CA', 'CB'), res='1')
    ':1@CA :1@CB'
    >>> join_mask('CA CB', res='1')
    ':1@CA :1@CB'
    >>> join_mask('CA CB', res=0)
    ':1@CA :1@CB'
    """
    from pytraj.compat import string_types

    if is_int(res):
        res = str(res + 1)
    else:
        res = res

    if isinstance(m, string_types):
        # 'CA CB' to ['CA', 'CB']
        m = m.split()
    elif not isinstance(m, (list, tuple)):
        raise ValueError("must be a list/tuple")

    return " ".join(':' + res + '@' + s for s in m)


def split_range(n_chunks, start, stop):
    '''split a given range to n_chunks

    Examples
    --------
    >>> split_range(3, 0, 10)
    [(0, 3), (3, 6), (6, 10)]
    '''
    list_of_tuple = []
    chunksize = (stop - start) // n_chunks
    for i in range(n_chunks):
        if i < n_chunks - 1:
            _stop = start + (i + 1) * chunksize
        else:
            _stop = stop
        list_of_tuple.append((start + i * chunksize, _stop))
    return list_of_tuple
