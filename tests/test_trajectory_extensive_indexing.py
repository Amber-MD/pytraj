import os
from time import time
import unittest
import pytraj as pt
import numpy as np
from pytraj.base import *
from pytraj.testing import aa_eq

refilename = "./data/Tc5b.nat.crd"
trajin = """
"""

ts = pt.iterload('data/Tc5b.x', 'data/Tc5b.top')

# create Trajectory to store Frame
FARRAY = Trajectory()
FRAMENUM = 999
FARRAY = ts[:FRAMENUM]


class TestTrajectory(unittest.TestCase):

    def test_len(self):
        N = 10
        farray = FARRAY[:N].copy()
        assert farray.n_frames == N
        old_xyz_5_10 = farray[5].xyz[:10]
        assert farray[:3].n_frames == 3
        assert farray[1:3].n_frames == 2
        assert farray[3:1].n_frames == 0
        assert farray[3:1:-1].n_frames == 2
        assert farray[-1:-3].n_frames == 0
        assert farray[-1:-3:-1].n_frames == 2
        aa_eq(farray[-1].xyz, farray[N - 1].xyz)

        # need to create a temp farray
        subfarray = farray[5:1:-1]
        aa_eq(subfarray[0].xyz, farray[5].xyz)
        aa_eq(old_xyz_5_10, farray[5].xyz[:10])

        f_last = farray[-3:-1][-1]

    def test_len_TrajectoryIterator(self):
        # create alias of `ts` (TrajectoryIterator instance  created above)
        farray = ts
        N = ts.n_frames
        assert farray.n_frames == N
        old_xyz_5_10 = farray[5].xyz[:10].copy()
        assert farray[:3].n_frames == 3
        assert farray[1:3].n_frames == 2
        assert farray[3:1].n_frames == 0
        assert farray[3:1:-1].n_frames == 2
        assert farray[-1:-3].n_frames == 0
        assert farray[-1:-3:-1].n_frames == 2
        # need to store xyz
        xyz = farray[-1].xyz.copy()
        aa_eq(xyz, farray[N - 1].xyz)

        # need to create a temp farray
        subfarray = farray[5:1:-1]
        aa_eq(subfarray[0].xyz, farray[5].xyz)
        aa_eq(old_xyz_5_10, farray[5].xyz[:10])

    def test_mask_indexing_0(self):
        # Trajectory
        traj = ts[:]
        assert traj["@CA"].shape == (traj.n_frames, traj.top("@CA").n_atoms, 3)
        assert traj[2:4]["@CA"].shape == (2, traj.top("@CA").n_atoms, 3)

    def test_mask_indexing(self):
        # Trajin_Single
        traj = ts
        assert traj["@CA"].shape == (traj.n_frames, traj.top("@CA").n_atoms, 3)
        assert traj[2:4]["@CA"].shape == (2, traj.top("@CA").n_atoms, 3)


if __name__ == "__main__":
    unittest.main()
