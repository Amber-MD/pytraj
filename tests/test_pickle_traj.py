from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.testing import eq, aa_eq, assert_equal_topology


class TestPickleTrajectoryIterator(unittest.TestCase):
    def test_trajiter(self):
        for frame_slice in [(0, 8, 2), (0, 10, 1)]:
            traj = pt.iterload("data/md1_prod.Tc5b.x", "data/Tc5b.top",
                    frame_slice=frame_slice)
            pt.io.to_pickle(traj, 'output/test0.pk')
            t0 = pt.io.read_pickle('output/test0.pk')

            aa_eq(traj.xyz, t0.xyz)
            assert_equal_topology(traj.top, t0.top, traj)


if __name__ == "__main__":
    unittest.main()
