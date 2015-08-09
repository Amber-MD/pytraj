from pytraj.testing import make_fake_traj, Timer, duplicate_traj
import pytraj.common_actions as pyca
import pytraj.io as io
from pytraj.misc import merge_trajs
from pytraj import Frame

traj = io.load_sample_data('tz2')[:]

for i in range(3):
    traj = merge_trajs(traj.copy(), traj.copy()) 

traj = duplicate_traj(traj, 300)
xyz = traj.xyz[:].copy()
f0 = traj[0].copy()
assert isinstance(f0, Frame) == True
xyz0 = f0.xyz.copy()

print (traj)

@Timer()
def test0(traj, f0):
    traj += f0

print ("traj")
test0(traj, f0)
print ("numpy")
test0(xyz, xyz0)
