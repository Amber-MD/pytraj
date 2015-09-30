import unittest
import pytraj as pt
import numpy as np
from array import array
import pytraj as pt
from pytraj.base import *
from pytraj.io import load
from pytraj.testing import aa_eq


class TestIndices(unittest.TestCase):
    def test_slice(self):

        traj1 = TrajectoryIterator(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")
        frame_indices = slice(9, 6, -1)

        traj0 = pt.load(filename="./data/md1_prod.Tc5b.x",
                        top=pt.load_topology("./data/Tc5b.top"),
                        frame_indices=frame_indices)

        assert traj0[0].same_coords_as(traj1[9]) == True
        assert traj0[1].same_coords_as(traj1[8]) == True
        assert traj0[2].same_coords_as(traj1[7]) == True

        assert traj0[0].rmsd(traj1[9]) < 1E-4

        rmsdlist = []
        ref = traj1[0].copy()
        for frame in traj1:
            rmsdlist.append(frame.rmsd(ref))

        nparr = np.array(rmsdlist)

        # make sure we don't suport other frame_indices
        traj2 = Trajectory()
        traj2 = pt.load(
            filename="./data/md1_prod.Tc5b.x",
            top=pt.load_topology("./data/Tc5b.top"),
            frame_indices=list(range(4)) + list(range(9, 5, -1)) + [4, ])
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
        top = pt.load_topology("./data/Tc5b.top")
        top2 = top
        del top

    def test_join_dummy(self):
        traj0 = TrajectoryIterator(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")[:]
        #traj0 += traj0
        traj0.join(traj0[:])

    def test_load_frame_indices_from_io(self):
        traj0 = pt.load(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top",
            frame_indices=(1, 3, 7))
        trajCA = pt.load(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top",
            frame_indices=(1, 3, 7),
            mask='@CA')
        trajreadonly = pt.iterload(
            filename="data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")
        trajCA_10frames = trajreadonly['@CA']

        assert isinstance(traj0, Trajectory)
        aa_eq(traj0[0].coords, trajreadonly[1].coords)
        aa_eq(traj0[1].coords, trajreadonly[3].coords)
        aa_eq(traj0[2].coords, trajreadonly[7].coords)

        # @CA
        aa_eq(trajCA[0].coords, trajCA_10frames[1].coords)
        aa_eq(trajCA[1].coords, trajCA_10frames[3].coords)
        aa_eq(trajCA[2].coords, trajCA_10frames[7].coords)


if __name__ == "__main__":
    unittest.main()
