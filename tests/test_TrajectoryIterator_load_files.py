from __future__ import print_function
import unittest
from pytraj import io
from pytraj.utils import eq, aa_eq

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import TrajectoryIterator
        traj = TrajectoryIterator("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        n_frame0 = 10
        assert traj.n_frames == n_frame0

        N = 3
        traj = TrajectoryIterator(["./data/md1_prod.Tc5b.x" for _ in range(N)],
                                   "./data/Tc5b.top")
        assert traj.n_frames == n_frame0 * N

        traj = TrajectoryIterator(["./data/md1_prod.Tc5b.x" for _ in range(N)],
                                   "./data/Tc5b.top", [(0, 5), (3, 8)])
        assert traj.n_frames == 10

if __name__ == "__main__":
    unittest.main()
