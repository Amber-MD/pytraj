from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        farray = traj[:]
        traj0_CA = traj[:]
        traj0_CA.strip("!@CA")

        # test TrajectoryIterator
        for idx, f0 in enumerate(traj(mask='@CA')):
            f1 = traj0_CA[idx]
            aa_eq(f0.xyz, f1.xyz)

        # test TrajectoryIterator with indices as mask
        indices = traj.top("@CA").indices
        for idx, f0 in enumerate(traj(mask=indices)):
            f1 = traj0_CA[idx]
            aa_eq(f0.xyz, f1.xyz)

        # test Trajectory
        for idx, f0 in enumerate(farray(mask='@CA')):
            f1 = traj0_CA[idx]
            aa_eq(f0.xyz, f1.xyz)

        # test Trajectory with indices
        for idx, f0 in enumerate(farray(mask=indices)):
            f1 = traj0_CA[idx]
            aa_eq(f0.xyz, f1.xyz)

        assert idx + 1 == traj0_CA.n_frames


if __name__ == "__main__":
    unittest.main()
