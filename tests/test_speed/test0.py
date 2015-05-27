from pytraj.testing import make_fake_traj, Timer, duplicate_traj
import pytraj.common_actions as pyca
import pytraj.io as io
from pytraj.misc import merge_trajs

#traj = make_fake_traj(10, 200000)
traj0 = io.load_sample_data('tz2') [:]
traj = traj0.copy()

traj.autoimage()
traj.rmsfit_to(0)

for i in range(5):
    traj = merge_trajs(traj.copy(), traj.copy()) 

traj = duplicate_traj(traj, 20)

@Timer()
def test0(mode='pytraj'):
    traj.rmsd(mode=mode, mask='!@H=')

print (traj)
print ("pytraj")
test0(mode="pytraj")
print ("cpptraj")
test0(mode="cpptraj")
