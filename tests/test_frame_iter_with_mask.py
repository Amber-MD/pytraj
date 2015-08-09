from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        print(farray)
        traj0_CA = traj[:]
        traj0_CA.strip_atoms("!@CA")
        print("after stripping")
        print(farray, traj0_CA, traj0_CA.n_frames)

        # test TrajectoryIterator
        for idx, f0 in enumerate(traj(mask='@CA')):
            print(idx, traj0_CA.n_frames)
            f1 = traj0_CA[idx]
            print(idx, f0, f1)
            assert_almost_equal(f0.coords, f1.coords)
            print(f0.rmsd(f1))

        # test TrajectoryIterator with indices as mask
        indices = traj.top("@CA").selected_indices()
        for idx, f0 in enumerate(traj(mask=indices)):
            f1 = traj0_CA[idx]
            print(idx, f0, f1)
            assert_almost_equal(f0.coords, f1.coords)
            print(f0.rmsd(f1))

        # test Trajectory
        for idx, f0 in enumerate(farray(mask='@CA')):
            f1 = traj0_CA[idx]
            print(idx, f0, f1)
            assert_almost_equal(f0.coords, f1.coords)
            print(f0.rmsd(f1))

        # test Trajectory with indices
        for idx, f0 in enumerate(farray(mask=indices)):
            f1 = traj0_CA[idx]
            print(idx, f0, f1)
            assert_almost_equal(f0.coords, f1.coords)
            print(f0.rmsd(f1))

        assert idx + 1 == traj0_CA.n_frames


if __name__ == "__main__":
    unittest.main()
