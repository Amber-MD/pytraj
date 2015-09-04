import unittest
import pytraj as pt
import numpy as np
from array import array
from pytraj.base import *
from pytraj.io import load
from pytraj.decorators import no_test
from pytraj.utils.check_and_assert import assert_almost_equal

class TestIndices(unittest.TestCase):
    #@no_test

    def test_0(self):

        traj1 = TrajectoryIterator(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")
        indices = slice(9, 6, -1)

        traj0 = pt.load(filename="./data/md1_prod.Tc5b.x",
                        top=Topology("./data/Tc5b.top"),
                        indices=indices)

        assert traj0[0].same_coords_as(traj1[9]) == True
        assert traj0[1].same_coords_as(traj1[8]) == True
        assert traj0[2].same_coords_as(traj1[7]) == True

        assert traj0[0].rmsd(traj1[9]) < 1E-4

        rmsdlist = []
        ref = traj1[0].copy()
        for frame in traj1:
            rmsdlist.append(frame.rmsd(ref))

        nparr = np.array(rmsdlist)

        # make sure we don't suport other indices
        traj2 = Trajectory()
        traj2 = pt.load(filename="./data/md1_prod.Tc5b.x",
                        top=Topology("./data/Tc5b.top"),
                        indices=list(range(4)) + list(range(9, 5, -1)) + [4, ])
        assert traj2[-1].coords == traj1[4].coords

    def test_array_assigment(self):
        traj1 = TrajectoryIterator(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")[:]

        # assign traj1[0]
        traj1[0] = traj1[1].copy()
        # make sure the assignment happed correctly
        assert traj1[0].same_coords_as(traj1[1]) == True

        traj1[0][10, 0] = 1000000.
        assert traj1[0][10, 0] == traj1[0, 10, 0] == 1000000.
        assert (traj1[0].same_coords_as(traj1[1])) == False
        assert traj1[0, 10, 0] != traj1[1, 10, 0]

    def test_1(self):
        traj0 = TrajectoryIterator(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")
        traj = TrajectoryIterator(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")[:]
        assert traj[0].coords == traj0[0].coords

        traj2 = TrajectoryIterator(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")[:][:10]
        assert traj2[0].coords == traj0[0].coords

        traj.join((traj[:], traj[0:100], traj[9:3:-1]))
        traj.join(traj[:])

        assert traj[0].coords != array('d', [0 for _ in range(traj[0].size)])
        assert traj[-1].coords != array('d', [0 for _ in range(traj[0].size)])

        for frame in traj:
            frame.zero_coords()

        assert traj[0].coords == array('d', [0 for _ in range(traj[0].size)])
        assert traj[-1].coords == array('d', [0 for _ in range(traj[0].size)])

    def test_del_top(self):
        # why here? lazy to make another file
        top = Topology("./data/Tc5b.top")
        top2 = top
        del top

    def test_join_dummy(self):
        traj0 = TrajectoryIterator(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")[:]
        #traj0 += traj0
        traj0.join(traj0[:])

    def test_load_indices_from_io(self):
        from pytraj import io as mdio
        traj0 = mdio.load(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top",
            indices=(1, 3, 7))
        trajreadonly = mdio.iterload(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")

        assert isinstance(traj0, Trajectory)
        assert_almost_equal(traj0[0].coords, trajreadonly[1].coords)
        assert_almost_equal(traj0[1].coords, trajreadonly[3].coords)
        assert_almost_equal(traj0[2].coords, trajreadonly[7].coords)


if __name__ == "__main__":
    unittest.main()
