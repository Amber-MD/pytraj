# Conclusion: traj.apply using `inplace` has similiar speed with numpy `inplace`
# without `inplace`, traj.apply is slightly faster than numpy
from pytraj.testing import make_fake_traj, Timer, duplicate_traj, aa_eq
import pytraj.common_actions as pyca
import pytraj.io as io
from pytraj.misc import merge_trajs

traj = io.load_sample_data('tz2')[:]

for i in range(5):
    traj = merge_trajs(traj.copy(), traj.copy()) 

traj = duplicate_traj(traj, 50)
xyz = traj.xyz[:].copy()

print (traj)

@Timer()
def test0(xyz, xyz0):
    xyz += xyz0

@Timer()
def test_apply(traj):
    traj.apply(lambda x : 2 * x + 1)

@Timer()
def test_apply_inplace(traj):
    def func(x):
        x *= 2.
        x += 1.
    traj.apply(func)

@Timer()
def test_apply_similiar_numpy_inplace(xyz):
    xyz *= 2.
    xyz += 1.

@Timer()
def test_apply_similiar_numpy(xyz):
    xyz = xyz * 2. + 1.

def run_0():
    f0 = traj[0].copy()
    print ("traj")
    test0(traj, f0)
    xyz = traj.xyz.copy()
    xyz0 = traj.xyz[0].copy()
    print ("numpy")
    test0(xyz, xyz0)

def run_apply_test():
    xyz = traj.xyz[:].copy()
    #traj.apply(lambda x : 2 * x + 1)
    #xyz = 2 * xyz + 1 
    #aa_eq(traj.xyz, xyz)

    print ("traj")
    test_apply(traj)
    print ("traj inplace")
    test_apply_inplace(traj)
    print ("numpy")
    test_apply_similiar_numpy(xyz)
    print ("numpy_inplace")
    test_apply_similiar_numpy_inplace(xyz)

run_apply_test()
