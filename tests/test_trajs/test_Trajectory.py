import unittest
import pytraj as pt
import numpy as np
from pytraj.base import *

traj = TrajectoryIterator()
traj.top = pt.load_topology("./data/Tc5b.top")
traj.load("./data/md1_prod.Tc5b.x")


class TestTrajectory(unittest.TestCase):
    def test_slice_basic(self):
        traj2 = Trajectory()
        traj2.top = pt.load_topology("./data/Tc5b.top")
        traj2.load("./data/md1_prod.Tc5b.x")
        traj2.load("./data/md1_prod.Tc5b.x")
        traj2.load("./data/md1_prod.Tc5b.x")
        traj2.load("./data/md1_prod.Tc5b.x")
        fsub = traj2[2:10]
        fsub[0][0] = 100.

    def test_indexing_0(self):
        traj2 = TrajectoryIterator()
        traj2.top = pt.load_topology("./data/Tc5b.top")
        traj2.load("./data/md1_prod.Tc5b.x")
        farray = traj2[[0, 9, 1]]
        assert farray.n_frames == 3
        assert traj2[0].atoms(0) == farray[0].atoms(0)
        assert traj2[9].atoms(0) == farray[1].atoms(0)
        assert traj2[1].atoms(0) == farray[2].atoms(0)

        arr = np.asarray(traj2[0].buffer1d[:])
        frame0 = traj2[0]
        arr0 = np.asarray(frame0.buffer1d[:])

        mat0 = np.asmatrix(arr0).reshape(304, 3)
        mat0[:, 0] = np.asmatrix(list(range(304))).reshape(304, 1)
        assert frame0[0, 0] == 0.
        assert frame0[1, 0] == 1.
        assert frame0[2, 0] == 2.

        # raise if size = 0
        traj3 = pt.Trajectory()
        assert traj3.n_frames == 0, 'empty Trajectory, n_frames must be 0'
        self.assertRaises(IndexError, lambda: traj3[0])
        self.assertRaises(IndexError, lambda: traj3.__setitem__(0, traj[3]))

    def test_indexing_1(self):
        traj2 = TrajectoryIterator()
        traj2.top = pt.load_topology("./data/Tc5b.top")
        traj2.load("./data/md1_prod.Tc5b.x")
        self.assertRaises(ValueError, lambda: traj2[10407])
        assert traj2[0] != traj2[9]

        assert traj2[-1].coords[0] == traj2[9].coords[0]

    def test_iter_basic(self):
        traj = TrajectoryIterator()
        traj.top = pt.load_topology("./data/Tc5b.top")
        traj.load("./data/md1_prod.Tc5b.x")
        for frame in traj:
            pass

    def test_trj_top(self):
        traj = TrajectoryIterator()
        assert traj.top.is_empty() == True
        traj.top = pt.load_topology("./data/Tc5b.top")
        assert traj.top.is_empty() == False
        traj.load("./data/md1_prod.Tc5b.x")

    def test_xyz(self):
        traj = pt.datafiles.load_tz2_ortho()
        t0 = traj[:]

        def set_xyz_not_c_contiguous():
            t0.xyz = np.asfortranarray(traj.xyz)

        def set_xyz_not_same_n_atoms():
            traj1 = pt.load_sample_data('ala3')
            t0.xyz = traj1.xyz

        # not c_contiguous
        self.assertRaises(TypeError, lambda: set_xyz_not_c_contiguous())
        self.assertRaises(ValueError, lambda: set_xyz_not_same_n_atoms())

if __name__ == "__main__":
    unittest.main()
