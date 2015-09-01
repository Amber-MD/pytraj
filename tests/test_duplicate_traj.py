from __future__ import print_function
import unittest
import pytraj as pt
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.testing import duplicate_traj
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        n_frames = traj.n_frames

        # TrajectoryIterator
        dtraj = duplicate_traj(traj, 2)
        assert dtraj.n_frames == 2 * traj.n_frames
        aa_eq(dtraj[:n_frames].xyz, traj.xyz)
        aa_eq(dtraj[n_frames:].xyz, traj.xyz)

        # Trajectory
        dtraj = duplicate_traj(traj.to_mutable_trajectory(), 2)
        assert dtraj.n_frames == 2 * traj.n_frames
        aa_eq(dtraj[:n_frames].xyz, traj.xyz)
        aa_eq(dtraj[n_frames:].xyz, traj.xyz)


if __name__ == "__main__":
    unittest.main()
