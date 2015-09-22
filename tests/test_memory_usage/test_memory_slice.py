from __future__ import print_function
from pytraj import io as mdio
from pytraj import Frame
from pytraj.compat import range
from memory_profiler import profile


@profile
def main():
    import pytraj.common_actions as pyca
    from pytraj import Trajectory
    print("create TrajectoryIterator")
    traj = mdio.iterload(
        "./data/nogit/tip3p/md.trj", "./data/nogit/tip3p/tc5bwat.top")
    print("create mutable Trajectory in memory")
    # try to slice
    fa = traj[:100]
    traj[:100]
    traj[:100]
    fa[:]
    fa[2:10:3]

    # try to make a copy
    fa.copy()

    # try to perform some analyses
    pyca.calc_rmsd(fa[:100])
    pyca.search_hbonds(fa[:100])
    pyca.calc_rmsd(traj[:100])
    pyca.calc_rmsd(traj(stop=101), top=traj.top, ref=traj[0])
    pyca.search_hbonds(traj[:100])
    pyca.search_hbonds(traj(stop=101), top=traj.top)

    # try to call `xyz` (copy)
    traj[:100].xyz
    traj[:100].xyz
    traj[:100].xyz
    fa[:100].xyz

    # try to join Trajectory
    fa += fa

    # try to append without copying
    fa.append(fa[0], copy=False)
    # try to append with copying
    fa.append(fa[0], copy=True)

    # try to load
    fa.load(fa[:2])

    # try to create frames
    Frame(fa.n_atoms)
    Frame(2 * fa.n_atoms)
    Frame(4 * fa.n_atoms)

    # try to create empty Trajectory
    fa1 = Trajectory()
    fa2 = Trajectory()
    # try to create n_frames=100 without allocate Frame
    fa3 = Trajectory(n_frames=100)
    fa3.top = traj.top.copy()
    # try to assign new Frame
    fa3[0] = traj[0]
    fa3.top = traj.top
    fa3[0].xyz
    # fa3.xyz # FIXME: DON'T NEED, segfault (because you have not created
    # other frames yet)
    print(fa3[0, 0])

    # try to append new frames
    for _ in range(100):
        fa1.append(Frame(fa.n_atoms))

    # try to iterate Trajectory
    for f in fa:
        f

    # try to iterate Trajectory by frame_iter
    for f in fa():
        f

    # try to iterate Trajectory by frame_iter with mask
    for f in fa(mask='!@H='):
        f

    # try to iterate Trajectory by frame_iter with mask
    for f in fa(mask=':WAT'):
        f

    # try to iterate whole TrajectoryIterator, stop at 1000 frames
    for f in traj(stop=1000):
        f

    # try to iterate whole TrajectoryIterator
    for f in traj:
        f

    # try to _allocate memory
    fa1._allocate(traj.n_frames / 10, traj.n_atoms)


if __name__ == "__main__":
    main()
