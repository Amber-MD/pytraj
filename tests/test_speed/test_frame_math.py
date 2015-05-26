from pytraj.testing import make_fake_traj, Timer
import pytraj.common_actions as pyca
import pytraj.io as io
from pytraj.misc import merge_trajs

# Conclusion: similiar speed to umpy
f0 = make_fake_traj(1, 500000)[0]
f1 = f0.copy()

xyz = f0.xyz.copy()
xyz1 = f1.xyz.copy()

@Timer()
def test0(f0, f1):
    f0 += f1

print ("frame")
test0(f0, f1)
del f0
del f1

print ("numpy")
test0(xyz, xyz1)
