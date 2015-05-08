import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.core import Box
from array import array as pyarray
from pytraj.decorators import test_if_having
from pytraj.testing import eq, aa_eq

class TestBox(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        frame0.box_crd()
        print(frame0.box)
        frame0.boxview[:] = pyarray('d', [0.0, 1.0, 2.0, 3.0, 4.0, 6.])
        print(frame0.box)
        print(frame0.box.type)
        frame0.set_nobox()
        print(frame0.box)

    def test_help(self):
        print (Box.all_box_types())

    def test_1(self):
        box = Box()
        box.set_trunc_oct()
        print(box)
        print(box.type)
        box.set_nobox()
        print(box)
        print(box.type)

        dummy = 100.
        box.data[0] = dummy
        assert box.data[0] == dummy
        assert box.tolist()[0] == dummy

    @test_if_having("numpy")
    def test_1(self):
        import numpy as np
        box = Box()
        arr0 = np.arange(6).astype(np.float64)
        box.data[:] = arr0
        print (box.tolist())

        for idx, x in enumerate(arr0):
            assert box.data[idx] == x

        # set Box for Frame
        f1 = Frame()
        f1.box = box
        eq(f1.box.tolist(), box.tolist())

        f2 = Frame()
        f2.box = box.to_ndarray()
        eq(f2.box.tolist(), box.tolist())

    def test_real_box(self):
        traj = mdio.load("./data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        trajiter = mdio.iterload("./data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        aa_eq(traj.top.box.tolist(), [35.2627796623, 41.8455476799, 36.168629529, 90.0, 90.0, 90.0], decimal=1)
        for frame in traj:
            assert frame.box.type == 'ortho'
            aa_eq(frame.box.tolist(), [35.2627796623, 41.8455476799, 36.168629529, 90.0, 90.0, 90.0], decimal=1)

        # test box_to_ndarray
        arr0 = traj.box_to_ndarray()
        arr1 = trajiter.box_to_ndarray()

        for b0, b1, frame in zip(arr0, arr1, trajiter):
            aa_eq(b0,  b1)
            aa_eq(b0,  frame.box.to_ndarray())

if __name__ == "__main__":
    unittest.main()
