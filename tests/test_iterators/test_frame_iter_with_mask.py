from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("data/Tc5b.x", "data/Tc5b.top")
        farray = traj[:]
        traj0_CA = traj[:]
        traj0_CA.strip("!@CA")

        # test TrajectoryIterator
        for idx, f0 in enumerate(traj(mask='@CA')):
            f1 = traj0_CA[idx]
            assert_almost_equal(f0.xyz, f1.xyz)

        # test TrajectoryIterator with indices as mask
        indices = traj.top("@CA").indices
        for idx, f0 in enumerate(traj(mask=indices)):
            f1 = traj0_CA[idx]
            assert_almost_equal(f0.xyz, f1.xyz)

        # test Trajectory
        for idx, f0 in enumerate(farray(mask='@CA')):
            f1 = traj0_CA[idx]
            assert_almost_equal(f0.xyz, f1.xyz)

        # test Trajectory with indices
        for idx, f0 in enumerate(farray(mask=indices)):
            f1 = traj0_CA[idx]
            assert_almost_equal(f0.xyz, f1.xyz)

        assert idx + 1 == traj0_CA.n_frames


if __name__ == "__main__":
    unittest.main()
