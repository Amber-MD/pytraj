import unittest
import pytraj as pt
import numpy as np
from pytraj import Trajectory
from pytraj.testing import aa_eq


class TestSlicingTrajectory(unittest.TestCase):

    def test_array_like(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        FA = traj[:]

        # slicing with list or array
        indices = [1, 2, 3]
        fa = traj[indices]
        fa2 = FA[indices]
        fa3 = traj[range(1, 4)]
        fa4 = FA[range(1, 4)]
        self.assertIsInstance(fa, Trajectory)
        # from TrajectoryIterator
        aa_eq(fa[0].xyz, traj[1].xyz)
        aa_eq(fa[1].xyz, traj[2].xyz)
        # from Trajectory
        aa_eq(fa2[1].xyz, traj[2].xyz)
        aa_eq(fa2[0].xyz, traj[1].xyz)

        # from "range"
        aa_eq(fa3[1].xyz, traj[2].xyz)
        aa_eq(fa3[0].xyz, traj[1].xyz)
        aa_eq(fa4[1].xyz, traj[2].xyz)
        aa_eq(fa4[0].xyz, traj[1].xyz)

    def test_atommask(self):
        # AtomMask
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        xyz = traj.xyz[:]
        atm = traj.top("@CA")
        indices = atm.indices

        aa_eq(fa[0, atm], fa[0][atm])
        aa_eq(traj[0, atm], fa[0][atm])
        aa_eq(traj[0, atm, 0], fa[0][atm, 0])
        aa_eq(traj[0, atm, 0], xyz[0][indices][0])


class Test1(unittest.TestCase):

    def test_0(self):
        # create Trajectory from Trajing_Single
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")[:]
        aa_eq(traj[3, 3], traj[3][3, :])
        frame1 = traj[1]
        aa_eq(frame1[0], traj[1][:, :][0])
        assert traj[0, 0, 0] == -16.492
        assert traj[:, :, 0][0, 0] == traj[0, 0, 0]
        f0 = traj[0]
        farr0 = traj[:2]

        fa = traj[2:4]

        # we don't support traj[:, idx] or traj[:, idx, idy] since this give wrong answer
        #  got ~0.0 value
        aa_eq(traj[:, 0, 0], np.asarray(traj[0][0]))

        for i in range(traj[0]._buffer2d.shape[0]):
            aa_eq(traj[:, :, 0][i], traj[0]._buffer2d[i])

        # slicing with mask
        atm = traj.top("@CA")
        traj[atm]
        traj[:, atm]


class TestSegmentationFault(unittest.TestCase):

    def test_0(self):
        # NOTE: no assert, just check for segfault
        traj = pt.load("./data/Tc5b.x", "./data/Tc5b.top")
        trajiter = pt.load("./data/Tc5b.x", "./data/Tc5b.top")
        atm = traj.top("@CA")
        f0 = traj[5]
        f0 = traj[0]
        f0.top = traj.top
        f0['@CA']
        traj[0, '@CA']

        f0 = traj[0, '@CA']
        f1 = traj['@CA'][0]
        assert pt.tools.rmsd(f0.xyz, f1.xyz) == 0.0


if __name__ == "__main__":
    unittest.main()
