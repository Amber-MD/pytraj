import pytraj as pt
from pytraj.testing import Timer
from memory_profiler import profile


@profile
def test():
    fname = 'data/nogit/remd/remd.x.000'
    topname = 'data/nogit/remd/myparm.parm7'
    traj = pt.trajectory.Trajectory(fname, top=topname)
    print(traj.xyz.dtype)

    for f in traj:
        pass

    pt.radgyr(traj)
    pt.molsurf(traj, '@CA')
    pt.rmsd(traj, mask='@CA')
    pt.multidihedral(traj)


@Timer()
def test_fastiter(traj):
    for f in traj:
        pass


fname = 'data/nogit/remd/remd.x.000'
topname = 'data/nogit/remd/myparm.parm7'
traj = pt.load(fname, topname)
print(traj)
test_fastiter(traj)
