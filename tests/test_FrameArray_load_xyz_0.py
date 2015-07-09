from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.compat import izip
from pytraj.decorators import test_if_having


class Test(unittest.TestCase):

    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print("creat Trajectory from 3D array")
        farray = Trajectory()
        farray.top = traj.top.copy()
        arr0 = traj.xyz
        print(arr0.shape)
        farray.append_xyz(arr0)
        for f0, f1 in izip(farray, traj):
            #print (f0, f1)
            assert_almost_equal(f0.coords, f1.coords)

        print("creat Trajectory from 1D array")
        farray2 = Trajectory()
        farray2.top = farray.top.copy()
        arr0 = traj.xyz.flatten()
        print(arr0.shape)
        farray2.append_xyz(arr0)

        for f0, f1 in izip(farray2, traj):
            #print (f0, f1)
            assert_almost_equal(f0.coords, f1.coords)

        print("creat Trajectory from 2D array")
        farray3 = Trajectory()
        farray3.top = farray.top.copy()
        arr0 = traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3)
        print(arr0.shape)
        farray3.append_xyz(arr0)

        for f0, f1 in izip(farray3, traj):
            #print (f0, f1)
            assert_almost_equal(f0.coords, f1.coords)

        print("creat Trajectory from 2D array of memview")
        farray4 = Trajectory()
        farray4.top = farray.top.copy()

        for frame in traj:
            f0 = Frame()
            f0.append_xyz(frame.buffer2d)
            farray4.append(f0)

        for f0, f1 in izip(farray4, traj):
            #print (f0, f1)
            assert_almost_equal(f0.coords, f1.coords)

        print("creat Trajectory from 2D array of ndarray")
        farray5 = Trajectory()
        farray5.top = farray.top.copy()

        farray5.append_ndarray(traj.xyz)
        print(farray5.size, traj.xyz.shape)

        i = 0
        for f0, f1 in izip(farray5, traj):
            #print (f0, f1)
            i += 1
            assert_almost_equal(f0.coords, f1.coords)
        assert i == traj.size

        print("creat Trajectory from frame_iter with mask")
        strip_top = farray5.top.strip_atoms("!@CA", copy=True)
        farray6 = Trajectory(traj(mask='@CA'), strip_top)
        farray5.strip_atoms('!@CA')
        print(farray5)

        for f0, f1 in izip(farray6, farray5):
            assert_almost_equal(f0.coords, f1.coords)

if __name__ == "__main__":
    unittest.main()
