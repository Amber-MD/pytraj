import unittest
import numpy as np
from array import array
import pytraj as pt
from pytraj.base import *
from pytraj.io import load
from pytraj.testing import aa_eq


class TestIndices(unittest.TestCase):

    def test_slice(self):

        traj1 = TrajectoryIterator(filename="data/Tc5b.x",
                                   top="./data/Tc5b.top")
        frame_indices = slice(9, 6, -1)

        traj0 = pt.load(filename="./data/Tc5b.x",
                        top=pt.load_topology("./data/Tc5b.top"),
                        frame_indices=frame_indices)

        aa_eq(traj0[0].xyz, traj1[9].xyz)
        aa_eq(traj0[1].xyz, traj1[8].xyz)
        aa_eq(traj0[2].xyz, traj1[7].xyz)

        assert traj0[0].rmsd(traj1[9]) < 1E-4

        rmsdlist = []
        ref = traj1[0].copy()
        for frame in traj1:
            rmsdlist.append(frame.rmsd(ref))

        nparr = np.array(rmsdlist)

        # make sure we don't suport other frame_indices
        traj2 = Trajectory()
        traj2 = pt.load(
            filename="./data/Tc5b.x",
            top=pt.load_topology("./data/Tc5b.top"),
            frame_indices=list(range(4)) + list(range(9, 5, -1)) + [4, ])
        aa_eq(traj2[-1].xyz, traj1[4].xyz)

    def test_del_top(self):
        # why here? lazy to make another file
        top = pt.load_topology("./data/Tc5b.top")
        top2 = top
        del top

    def test_load_frame_indices_from_io(self):
        traj0 = pt.load(filename="data/Tc5b.x",
                        top="./data/Tc5b.top",
                        frame_indices=(1, 3, 7))
        trajCA = pt.load(filename="data/Tc5b.x",
                         top="./data/Tc5b.top",
                         frame_indices=(1, 3, 7),
                         mask='@CA')
        trajreadonly = pt.iterload(filename="data/Tc5b.x",
                                   top="./data/Tc5b.top")
        trajCA_10frames = trajreadonly['@CA']

        assert isinstance(traj0, Trajectory)
        aa_eq(traj0[0].xyz, trajreadonly[1].xyz)
        aa_eq(traj0[1].xyz, trajreadonly[3].xyz)
        aa_eq(traj0[2].xyz, trajreadonly[7].xyz)

        # @CA
        aa_eq(trajCA[0].xyz, trajCA_10frames[1].xyz)
        aa_eq(trajCA[1].xyz, trajCA_10frames[3].xyz)
        aa_eq(trajCA[2].xyz, trajCA_10frames[7].xyz)

    def test_load_mask(self):
        traj = pt.iterload(filename="data/Tc5b.x",
                           top="./data/Tc5b.top")
        t0 = pt.load(traj.filename, traj.top.filename, mask='@CA')
        aa_eq(traj['@CA'].xyz, t0.xyz)


if __name__ == "__main__":
    unittest.main()
