from pytraj import io
from pytraj.utils import Timer
import pytraj.common_actions as pyca
from pytraj._shared_methods import _frame_iter_master
import numpy as np
from numpy.testing import assert_almost_equal as aa_equal

traj = io.iterload("remd.x.000", "myparm.top")

#trajlist = [traj, traj]
# 200 Mb each
trajlist = [traj, ]

print("not using frame_iter")


@Timer()
def test_time():
    return pyca.calc_pairwise_rmsd(trajlist, '@CA', traj.top)

dslist = test_time()

print("not using frame_iter")
#traj_iter_list = [traj(mask='@CA'), traj(mask='@CA')]
traj_iter_list = [traj(mask='@CA'), ]

# ned to get stripped-atom top
new_top = traj.top.strip_atoms('!@CA', copy=True)


@Timer()
def test_time_iter():
    # we already specify mask in frame_iter
    # and need to specify top too
    return pyca.calc_pairwise_rmsd(traj_iter_list,
                                   top=new_top)

dslist_iter = test_time_iter()

#print (dslist.size)
#print (dslist_iter.size)
#arr = dslist[0].to_ndarray().flatten()
#arr_iter = dslist_iter[0].to_ndarray().flatten()
#aa_equal(arr, arr_iter)
