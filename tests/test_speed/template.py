from pytraj.testing import make_fake_traj, Timer, duplicate_traj
import pytraj.common_actions as pyca
import pytraj.io as io
from pytraj.misc import merge_trajs

traj = io.load_sample_data('tz2')[:]

for i in range(5):
    traj = merge_trajs(traj.copy(), traj.copy()) 

traj = duplicate_traj(traj, 20)
print (traj)

@Timer()
def test0():
    pass
