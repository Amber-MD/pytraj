import os
from time import time
import unittest
import numpy as np
from pytraj.base import *
from pytraj.Timer import Timer
from load_traj import load
from pytraj.decorators import no_test

ts = TrajReadOnly()
datadir = "./data/"
topname = datadir + "Tc5b.top"
refilename = "./data/Tc5b.nat.crd"
mdx = "./data/md1_prod.Tc5b.x"
ts = TrajReadOnly()
print(ts.top.tag)

top = Topology(topname)
trajin = """
"""

ts.load(mdx, top)
ts.prepare_for_read(True)
frame = Frame()
frame.set_frame_v(top, ts.has_vel(), ts.n_repdims)
frame2 = Frame(frame)

# create FrameArray to store Frame
FARRAY = FrameArray()
#FARRAY.get_frames(ts, update_top=True)
FRAMENUM=999
FARRAY = ts[:FRAMENUM]

class TestFrameArray(unittest.TestCase):
    def test_len(self):
        N = 10
        farray = FARRAY[:N].copy()
        assert farray.size == N
        old_coords_5_10 = farray[5].coords[:10]
        assert farray[:3].size == 3
        assert farray[1:3].size == 2
        assert farray[3:1].size == 0
        assert farray[3:1:-1].size == 2
        assert farray[-1:-3].size == 0
        assert farray[-1:-3:-1].size == 2
        assert farray[-1].same_coords_as(farray[N-1]) == True

        #assert farray[5:1:-1][0].same_coords_as(farray[5]) == True
        # segment fault if using below expression
        #print farray[5:1:-1][0].coords[:10]

        # need to create a temp farray
        subfarray = farray[5:1:-1]
        print(subfarray)
        assert subfarray[0].same_coords_as(farray[5]) == True
        assert old_coords_5_10 == farray[5].coords[:10]
        print(subfarray[0].coords[:10])
        print(farray[5].coords[:10])

        f_last = farray[-3:-1][-1]
        print("***********XXXXXXXXXXXXX*")
        print(f_last)
        #print f_last.coords[:10]
        #print farray[-1].coords[:10]
        #print farray[-2].coords[:10]
        #print farray[-3].coords[:10]
        #assert f_last.same_coords_as(farray[-2]) == True

    def test_len_TrajReadOnly(self):
        # create alias of `ts` (TrajReadOnly instance  created above)
        farray = ts
        N = ts.size
        assert farray.size == N
        old_coords_5_10 = farray[5].coords[:10]
        assert farray[:3].size == 3
        assert farray[1:3].size == 2
        assert farray[3:1].size == 0
        assert farray[3:1:-1].size == 2
        assert farray[-1:-3].size == 0
        assert farray[-1:-3:-1].size == 2
        assert farray[-1].same_coords_as(farray[N-1]) == True

        #assert farray[5:1:-1][0].same_coords_as(farray[5]) == True
        # segment fault if using below expression
        #print farray[5:1:-1][0].coords[:10]

        # need to create a temp farray
        subfarray = farray[5:1:-1]
        print(subfarray)
        assert subfarray[0].same_coords_as(farray[5]) == True
        assert old_coords_5_10 == farray[5].coords[:10]
        print(subfarray[0].coords[:10])
        print(farray[5].coords[:10])

    def test_mask_indexing_0(self):
        # FrameArray
        traj = ts[:]
        print(type(traj["@CA"]))
        print(traj["@CA"].shape)
        assert traj["@CA"].shape == (traj.size, traj.top("@CA").n_selected, 3)
        print(traj[2:4]["@CA"])
        assert traj[2:4]["@CA"].shape == (2, traj.top("@CA").n_selected, 3)

    def test_mask_indexing(self):
        # Trajin_Single
        traj = ts
        print(type(traj["@CA"]))
        print(traj["@CA"].shape)
        assert traj["@CA"].shape == (traj.size, traj.top("@CA").n_selected, 3)
        print(traj[2:4]["@CA"])
        assert traj[2:4]["@CA"].shape == (2, traj.top("@CA").n_selected, 3)

if __name__ == "__main__":
    unittest.main()
