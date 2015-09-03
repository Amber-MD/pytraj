import unittest
import pytraj as pt
import numpy as np
from pytraj.base import *
from pytraj.decorators import no_test

TRAJ = TrajectoryIterator()
TRAJ.top = Topology("./data/Tc5b.top")
TRAJ.load("./data/md1_prod.Tc5b.x")


class TestTrajectory(unittest.TestCase):
    def test_slice_basic(self):
        TRAJ2 = Trajectory()
        TRAJ2.top = Topology("./data/Tc5b.top")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        fsub = TRAJ2[2:10]
        fsub[0][0] = 100.

    def test_indexing_0(self):
        TRAJ2 = TrajectoryIterator()
        TRAJ2.top = Topology("./data/Tc5b.top")
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
        TRAJ2.top = Topology("./data/Tc5b.top")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        self.assertRaises(ValueError, lambda: TRAJ2[10407])
        assert TRAJ2[0] != TRAJ2[9]

        assert TRAJ2[-1].coords[0] == TRAJ2[9].coords[0]

    def test_iter_basic(self):
        TRAJ = TrajectoryIterator()
        TRAJ.top = Topology("./data/Tc5b.top")
        TRAJ.load("./data/md1_prod.Tc5b.x")
        for frame in TRAJ:
            pass

    def test_iter(self):
        farray = Trajectory()
        farray.top = TRAJ.top
        i = 0
        for frame in TRAJ:
            i += 1
            frame.strip_atoms(top=TRAJ.top.copy(), mask="!@CA")
            farray.append(frame.copy())
        assert i == TRAJ.n_frames == TRAJ.n_frames
        assert frame.size == TRAJ.top.n_residues * 3
        farray.top.strip_atoms("!@CA")
        assert farray.top.n_atoms == TRAJ.top.n_residues
        farray.top.summary()
        assert farray.n_frames == TRAJ.n_frames
        arr = np.zeros(farray.n_frames)
        cpptraj_rmsd = np.loadtxt(
            "./data/rmsd_to_firstFrame_CA_allres.Tc5b.dat",
            skiprows=1).transpose()[1]

        # caculate rmsd to 1st frame
        for i in range(farray.n_frames):
            arr[i] = farray[0].rmsd(farray[i])
        np.testing.assert_almost_equal(arr[:10], cpptraj_rmsd[:10], decimal=3)

    def test_trj_top(self):
        traj = TrajectoryIterator()
        assert traj.top.is_empty() == True
        traj.top = Topology("./data/Tc5b.top")
        assert traj.top.is_empty() == False
        traj.load("./data/md1_prod.Tc5b.x")

    def test_1(self):
        traj = TrajectoryIterator()

        traj.load("./data/md1_prod.Tc5b.x", Topology("./data/Tc5b.top"))

        traj.top.summary()
        assert traj.top.n_atoms == 304


if __name__ == "__main__":
    unittest.main()
