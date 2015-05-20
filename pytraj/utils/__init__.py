""""""
from __future__ import absolute_import
from .check_and_assert import assert_almost_equal, file_exist, is_generator
from .check_and_assert import assert_almost_equal as aa_eq
from .check_and_assert import eq, a_isinstance, eq_coords
from .check_and_assert import _import, _import_numpy, is_int, require
from .check_and_assert import has_, is_array
from .check_and_assert import ensure_not_none_or_string
from .list_misc import flatten_list
from .Timer import Timer
from .context import goto_temp_folder
from ..compat import range

# add amberhome
from .amber_test import amberhome, cpptraj_test_dir, has_amberhome

def duplicate_traj(orig_traj, n_times):
    # always make copy
    if 'Iter' in orig_traj.__class__.__name__:
        traj = orig_traj[:]
    else:
        traj = orig_traj

    fa = traj.copy()
    for _ in range(n_times-1):
        fa.join(traj.copy())
    return fa
