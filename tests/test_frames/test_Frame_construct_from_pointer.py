from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        for xyz0 in traj.xyz:
            frame = pt.Frame(traj.n_atoms, xyz0, _as_ptr=True)
            aa_eq(frame.xyz, xyz0)

    def test_1(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        trajectory_t0 = pt.trajectory.Trajectory(traj)

        # __iter__
        for f in trajectory_t0:
            pass

        f.xyz[0, 0] = 10.
        assert f.xyz[0, 0] == 10.
        assert trajectory_t0.xyz[-1, 0, 0] == 10.

        # __getitem__
        # make a new copy
        trajectory_t0 = pt.trajectory.Trajectory(traj)
        f0 = trajectory_t0[0]
        f0.xyz[0, 0] = 200.
        assert trajectory_t0.xyz[0, 0, 0] == 200.

        # translate
        # make a new copy
        trajectory_t0 = pt.trajectory.Trajectory(traj)
        t0 = traj[:]
        pt.translate(trajectory_t0, 'x 1.2')
        pt.translate(t0, 'x 1.2')
        aa_eq(trajectory_t0.xyz, t0.xyz)

        try:
            import mdtraj as md
            pt.rmsd(t0, ref=0)
            pt.rmsd(trajectory_t0, ref=0)
        except ImportError:
            pass

        # test rmsfit
        trajectory_t0.rmsfit(ref=0)


if __name__ == "__main__":
    unittest.main()
