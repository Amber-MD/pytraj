import pytraj as pt
from memory_profiler import profile

@profile
def test():
    fname = 'data/nogit/remd/remd.x.000'
    topname = 'data/nogit/remd/myparm.parm7'
    traj = pt.api.Trajectory(fname, top=topname)
    print(traj.xyz.dtype)

    for f in traj:
        pass

    pt.radgyr(traj)
    pt.molsurf(traj, '@CA')
    pt.rmsd(traj, mask='@CA')
    pt.multidihedral(traj)

@profile
def test_fastiter():
    fname = 'data/nogit/remd/remd.x.000'
    topname = 'data/nogit/remd/myparm.parm7'
    t0 = pt.load(fname, topname)
    traj = pt.api.Trajectory(t0)
    for f in traj:
        print(f.box)

#test()
test_fastiter()
