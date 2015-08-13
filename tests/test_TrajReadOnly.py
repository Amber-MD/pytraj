import unittest
import numpy as np
from pytraj.base import *
from pytraj.decorators import no_test

TRAJ = TrajectoryIterator()
TRAJ.top = Topology("./data/Tc5b.top")
TRAJ.load("./data/md1_prod.Tc5b.x")
print("TRAJ.size", TRAJ.size)


class TestTrajectory(unittest.TestCase):
    def test_slice_basic(self):
        print("test_slice_basic")
        TRAJ2 = Trajectory()
        print(TRAJ2)
        #TRAJ2.debug = True
        TRAJ2.top = Topology("./data/Tc5b.top")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        print("TRAJ.size", TRAJ2.size)
        print(TRAJ2[3])
        print(TRAJ2[1])
        print(TRAJ2[2].coords[:10])
        print(TRAJ2[2].coords[:10])
        print(TRAJ2[2].coords[:10])
        print(TRAJ2[2].coords[:10])
        print(TRAJ2[2].coords[:10])
        print(TRAJ2[7])
        #fsub = TRAJ2[2:10].copy()
        fsub = TRAJ2[2:10]
        print(fsub[0].coords[0])
        print(fsub[0].coords[0])
        fsub[0][0] = 100.
        print(fsub[0][0])
        print(TRAJ2[2].coords[0])
        # print TRAJ2[2:10].copy()[0].coords[:10]
        # print TRAJ2[::]
        #farray = TRAJ2[::]
        # print farray
        # print "REALLY REALLY "
        # print farray.size
        # print farray[:4][0].n_atoms
        # print farray[:200][0].n_atoms
        # print farray[:50][0].n_atoms
        #fsub = farray[:200]
        # print fsub[0].n_atoms

        print(TRAJ2[::][0].n_atoms)

    #@no_test
    def test_slice(self):
        print("test_slice")
        TRAJ2 = TrajectoryIterator()
        #TRAJ2.debug = True
        TRAJ2.top = Topology("./data/Tc5b.top")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        print(TRAJ2.size)
        farray0 = TRAJ2[0:3:1]

        print("after slicing")
        print("n_atoms = ", TRAJ2[9:1000:50][0].n_atoms)
        print(farray0.size)
        print(farray0[0].n_atoms)

        print("test dummy slicing")
        print(TRAJ2[0:3])
        farray1 = TRAJ2[0:3]
        print(farray1[0].n_atoms)
        print(farray1)
        print("end test dummy slicing")

        print("test dummy slicing 2")
        print(TRAJ2[:3])
        print(TRAJ2[:3:1])
        print(TRAJ2[:3:2])
        print(TRAJ2[:3:-1])
        print(TRAJ2[::1000])
        print(TRAJ2[1::10])

        print(TRAJ2[100:0:-1])
        print(TRAJ2[::1])
        print(TRAJ2[::])

        farray2 = Trajectory()
        indices = list(range(TRAJ2.size))
        farray2.get_frames(TRAJ2, update_top=True)
        farray3 = Trajectory()
        farray3.get_frames(farray2, update_top=True)

    #@no_test
    def test_indexing_0(self):
        print("test_indexing_0")
        TRAJ2 = TrajectoryIterator()
        TRAJ2.top = Topology("./data/Tc5b.top")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        farray = Trajectory()
        farray.get_frames(TRAJ2, indices=(0, 9, 1), update_top=True)
        print(farray.top)
        print(farray.size)
        print(TRAJ2.size)
        assert farray.size == 3
        print("************XDDFDFDFDFD")
        print(farray.size)
        assert TRAJ2[0].atoms(0) == farray[0].atoms(0)
        assert TRAJ2[9].atoms(0) == farray[1].atoms(0)
        assert TRAJ2[1].atoms(0) == farray[2].atoms(0)

        arr = np.asarray(TRAJ2[0].buffer1d[:])
        print("len")
        print(len(arr))
        print(TRAJ2[0].coords[:10])
        print("arr[:10]: ", arr[:10])
        print("test buffer1d view")
        frame0 = TRAJ2[0]
        arr0 = np.asarray(frame0.buffer1d[:])
        print("arr0[:10]: ", arr0[:10])
        print(TRAJ2[0].coords[:10])
        print(frame0.coords[:10])

        mat0 = np.asmatrix(arr0).reshape(304, 3)
        print(mat0.shape)
        mat0[:, 0] = np.asmatrix(list(range(304))).reshape(304, 1)
        assert frame0[0, 0] == 0.
        assert frame0[1, 0] == 1.
        assert frame0[2, 0] == 2.

    #@no_test
    def test_indexing_1(self):
        print("test_indexing_1")
        TRAJ2 = TrajectoryIterator()
        TRAJ2.top = Topology("./data/Tc5b.top")
        # TRAJ2.top.strip_atoms("!@CA")
        TRAJ2.load("./data/md1_prod.Tc5b.x")
        print(TRAJ2.size)
        print(TRAJ2.top.n_atoms)
        print(TRAJ2[0].coords[:10])
        print(TRAJ2[7].coords[:10])
        print("out of index testing")
        self.assertRaises(ValueError, lambda: TRAJ2[10407])
        assert TRAJ2[0] != TRAJ2[9]

        print("test negative indexing")
        print(TRAJ2[-1].coords[0])
        print(TRAJ2[9].coords[0])
        assert TRAJ2[-1].coords[0] == TRAJ2[9].coords[0]

    #@no_test
    def test_iter_basic(self):
        TRAJ = TrajectoryIterator()
        TRAJ.top = Topology("./data/Tc5b.top")
        TRAJ.load("./data/md1_prod.Tc5b.x")
        print("test_iter_basic")
        for frame in TRAJ:
            pass

    #@no_test
    def test_iter(self):
        farray = Trajectory()
        farray.top = TRAJ.top
        print("test_info")
        i = 0
        for frame in TRAJ:
            i += 1
            frame.strip_atoms(top=TRAJ.top.copy(), mask="!@CA")
            farray.append(frame.copy())
        assert i == TRAJ.size == TRAJ.n_frames
        assert frame.size == TRAJ.top.n_res * 3
        farray.top.strip_atoms("!@CA")
        print("farray.top.n_atoms= ", farray.top.n_atoms)
        assert farray.top.n_atoms == TRAJ.top.n_res
        farray.top.summary()
        assert farray.size == TRAJ.n_frames
        print("rmsd to first = ", farray[0].rmsd(farray[1]))
        arr = np.zeros(farray.size)
        cpptraj_rmsd = np.loadtxt(
            "./data/rmsd_to_firstFrame_CA_allres.Tc5b.dat",
            skiprows=1).transpose()[1]
        print(cpptraj_rmsd[:10])

        # caculate rmsd to 1st frame
        for i in range(farray.size):
            arr[i] = farray[0].rmsd(farray[i])
        print(arr[:10])
        np.testing.assert_almost_equal(arr[:10], cpptraj_rmsd[:10], decimal=3)
        print("Kool, reproduce cpptraj output")

    #@no_test
    def test_trj_top(self):
        traj = TrajectoryIterator()
        print(traj.top.is_empty())
        assert traj.top.is_empty() == True
        traj.top = Topology("./data/Tc5b.top")
        # traj.top.summary()
        assert traj.top.is_empty() == False
        traj.load("./data/md1_prod.Tc5b.x")
        #traj.load("./data/md1_prod.Tc5b.x", Topology("./data/Tc5b.top"))
        print(traj.size)

    #@no_test
    def test_1(self):
        traj = TrajectoryIterator()

        traj.load("./data/md1_prod.Tc5b.x", Topology("./data/Tc5b.top"))

        traj.top.summary()
        assert traj.top.n_atoms == 304
        print(traj.top.n_atoms)


if __name__ == "__main__":
    unittest.main()
