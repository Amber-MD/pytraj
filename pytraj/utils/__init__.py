""""""
from __future__ import absolute_import
from .check_and_assert import assert_almost_equal, file_exist, is_generator
from .check_and_assert import assert_almost_equal as aa_eq
from .check_and_assert import eq, a_isinstance, eq_coords
from .check_and_assert import _import, _import_numpy, is_int, require, _import_pandas
from .check_and_assert import has_, is_array
from .check_and_assert import ensure_not_none_or_string
from .list_misc import flatten_list
from .Timer import Timer
from .context import goto_temp_folder
from ..compat import range
from .import convert

# add amberhome
from .amber_test import amberhome, cpptraj_test_dir, has_amberhome


def duplicate_traj(orig_traj, n_times):
    # always make copy
    traj = orig_traj.copy()
    for _ in range(n_times - 1):
        if 'Iter' in orig_traj.__class__.__name__:
            # TrajectoryIterator
            traj.load(orig_traj.filelist)
        else:
            # Trajectory
            traj.join(traj.copy())
    return traj

def join_mask(m, res=None):
    """ join_mask(('CA', 'CB'), res='1') return ':1@CA :1@CB'
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
