import os
from time import time
import unittest
import pytraj as pt
import numpy as np
from pytraj.base import *
from pytraj.decorators import no_test

ts = TrajectoryIterator()
datadir = "./data/"
topname = datadir + "Tc5b.top"
refilename = "./data/Tc5b.nat.crd"
mdx = "./data/md1_prod.Tc5b.x"
ts = TrajectoryIterator()

top = pt.load_topology(topname)
trajin = """
"""

ts.load(mdx, top)
frame = Frame()
frame.set_frame_v(top)
frame2 = Frame(frame)

# create Trajectory to store Frame
FARRAY = Trajectory()
#FARRAY.get_frames(ts, update_top=True)
FRAMENUM = 999
FARRAY = ts[:FRAMENUM]


class TestTrajectory(unittest.TestCase):
    def test_len(self):
        N = 10
        farray = FARRAY[:N].copy()
        assert farray.n_frames == N
        old_coords_5_10 = farray[5].coords[:10]
        assert farray[:3].n_frames == 3
        assert farray[1:3].n_frames == 2
        assert farray[3:1].n_frames == 0
        assert farray[3:1:-1].n_frames == 2
        assert farray[-1:-3].n_frames == 0
        assert farray[-1:-3:-1].n_frames == 2
        assert farray[-1].same_coords_as(farray[N - 1]) == True

        #assert farray[5:1:-1][0].same_coords_as(farray[5]) == True
        # segment fault if using below expression

        # need to create a temp farray
        subfarray = farray[5:1:-1]
        assert subfarray[0].same_coords_as(farray[5]) == True
        assert old_coords_5_10 == farray[5].coords[:10]

        f_last = farray[-3:-1][-1]
        #assert f_last.same_coords_as(farray[-2]) == True

    def test_len_TrajectoryIterator(self):
        # create alias of `ts` (TrajectoryIterator instance  created above)
        farray = ts
        N = ts.n_frames
        assert farray.n_frames == N
        old_coords_5_10 = farray[5].coords[:10]
        assert farray[:3].n_frames == 3
        assert farray[1:3].n_frames == 2
        assert farray[3:1].n_frames == 0
        assert farray[3:1:-1].n_frames == 2
        assert farray[-1:-3].n_frames == 0
        assert farray[-1:-3:-1].n_frames == 2
        assert farray[-1].same_coords_as(farray[N - 1]) == True

        #assert farray[5:1:-1][0].same_coords_as(farray[5]) == True
        # segment fault if using below expression

        # need to create a temp farray
        subfarray = farray[5:1:-1]
        assert subfarray[0].same_coords_as(farray[5]) == True
        assert old_coords_5_10 == farray[5].coords[:10]

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
