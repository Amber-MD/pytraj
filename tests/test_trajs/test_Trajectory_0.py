import unittest
import pytraj as pt
import numpy as np
from pytraj.base import *

TRAJ = TrajectoryIterator()
TRAJ.top = pt.load_topology("./data/Tc5b.top")
TRAJ.load("./data/md1_prod.Tc5b.x")


class TestTrajectory(unittest.TestCase):
    def test_slice_basic(self):
        TRAJ2 = Trajectory()
        TRAJ2.top = pt.load_topology("./data/Tc5b.top")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        fsub = TRAJ2[2:10]
        fsub[0][0] = 100.

    def test_indexing_0(self):
        TRAJ2 = TrajectoryIterator()
        TRAJ2.top = pt.load_topology("./data/Tc5b.top")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        farray = TRAJ2[[0, 9, 1]]
        assert farray.n_frames == 3
        assert TRAJ2[0].atoms(0) == farray[0].atoms(0)
        assert TRAJ2[9].atoms(0) == farray[1].atoms(0)
        assert TRAJ2[1].atoms(0) == farray[2].atoms(0)

        arr = np.asarray(TRAJ2[0].buffer1d[:])
        frame0 = TRAJ2[0]
        arr0 = np.asarray(frame0.buffer1d[:])

        mat0 = np.asmatrix(arr0).reshape(304, 3)
        mat0[:, 0] = np.asmatrix(list(range(304))).reshape(304, 1)
        assert frame0[0, 0] == 0.
        assert frame0[1, 0] == 1.
        assert frame0[2, 0] == 2.

    def test_indexing_1(self):
        TRAJ2 = TrajectoryIterator()
        TRAJ2.top = pt.load_topology("./data/Tc5b.top")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        self.assertRaises(ValueError, lambda: TRAJ2[10407])
        assert TRAJ2[0] != TRAJ2[9]

        assert TRAJ2[-1].coords[0] == TRAJ2[9].coords[0]

    def test_iter_basic(self):
        TRAJ = TrajectoryIterator()
        TRAJ.top = pt.load_topology("./data/Tc5b.top")
        TRAJ.load("./data/md1_prod.Tc5b.x")
        for frame in TRAJ:
            pass

    def test_trj_top(self):
        traj = TrajectoryIterator()
        assert traj.top.is_empty() == True
        traj.top = pt.load_topology("./data/Tc5b.top")
        assert traj.top.is_empty() == False
        traj.load("./data/md1_prod.Tc5b.x")

    def test_1(self):
        traj = TrajectoryIterator()

        traj.load("./data/md1_prod.Tc5b.x", pt.load_topology(
            "./data/Tc5b.top"))
        assert traj.top.n_atoms == 304


if __name__ == "__main__":
    unittest.main()
