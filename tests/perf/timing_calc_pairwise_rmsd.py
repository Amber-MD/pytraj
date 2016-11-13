from pytraj import io
from pytraj.utils import Timer


traj = io.iterload("remd.x.000", "myparm.top")

# 200 Mb each
trajlist = [traj, ]

print("not using frame_iter")


@Timer()
def test_time():
    return pt.calc_pairwise_rmsd(trajlist, '@CA', traj.top)


dslist = test_time()

print("not using frame_iter")
traj_iter_list = [traj(mask='@CA'), ]

# ned to get stripped-atom top
new_top = traj.top.strip('!@CA', copy=True)


@Timer()
def test_time_iter():
    # we already specify mask in frame_iter
    # and need to specify top too
    return pt.calc_pairwise_rmsd(traj_iter_list, top=new_top)


dslist_iter = test_time_iter()
