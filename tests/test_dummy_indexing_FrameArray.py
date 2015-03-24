import os
import unittest
import numpy as np
from pytraj.base import *
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
frame = Frame()
frame.set_frame_v(top)
frame2 = Frame(frame)

# create FrameArray to store Frame
FARRAY = FrameArray()
#FARRAY.get_frames(ts, update_top=True)
FRAMENUM=1000
FARRAY = ts[:FRAMENUM]

class TestFrameArray(unittest.TestCase):
    def test_dummy(self):
        farray = FARRAY.copy()
        print(farray[:])
        print(farray[1:4])
        print(farray[1:4])
        print(farray[1:4:2])
        print("farray[::-1]", farray[::-1])
        print("farray[::]", farray[::])
        print(farray[-1])
        print(farray[-1:-3])
        print(farray[-3:-1])
        print(farray[-3:-1])
        print(farray[-1:-10:2])
        print(farray[-1:-10:-2])
        print(farray[-1:-10:-2])
        newfarray = farray[:5] + farray[10:20]
        newfarray[0][0] = 100000.
        assert newfarray[0][0, 0] == 100000.
        print(farray[0][0])
        print("farray[::-1]", farray[::-1])
        farray2 = farray[::-1] + farray[:]
        print(farray2)
        print(farray2[0].n_atoms)

if __name__ == "__main__":
    unittest.main()
