from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.six_2 import izip

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print ("creat FrameArray from 3D array")
        farray = FrameArray()
        farray.top = traj.top.copy()
        arr0 = traj.xyz
        print (arr0.shape)
        farray.load_xyz(arr0)
        for f0, f1 in izip(farray, traj):
            print (f0, f1)
            assert_almost_equal(f0.coords, f1.coords)

        print ("creat FrameArray from 1D array")
        farray2 = FrameArray()
        farray2.top = farray.top.copy()
        arr0 = traj.xyz.flatten()
        print (arr0.shape)
        farray2.load_xyz(arr0)

        for f0, f1 in izip(farray2, traj):
            print (f0, f1)
            assert_almost_equal(f0.coords, f1.coords)

        print ("creat FrameArray from 2D array")
        farray3 = FrameArray()
        farray3.top = farray.top.copy()
        arr0 = traj.xyz.reshape(traj.n_frames, traj.n_atoms*3)
        print (arr0.shape)
        farray3.load_xyz(arr0)

        for f0, f1 in izip(farray3, traj):
            print (f0, f1)
            assert_almost_equal(f0.coords, f1.coords)

if __name__ == "__main__":
    unittest.main()
