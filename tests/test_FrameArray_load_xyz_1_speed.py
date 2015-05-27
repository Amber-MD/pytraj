from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.compat import izip
from pytraj.decorators import test_if_having
from pytraj.utils import Timer

class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        from numpy.testing import assert_almost_equal
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        @Timer()
        def _f_ndarray(farray, xyz):
            farray.append_ndarray(xyz)

        @Timer()
        def _f_dontknow(farray, xyz):
            farray.append_xyz(xyz)

        xyz = traj.xyz
        print (xyz.shape)

        f0 = Trajectory()
        f0.top = traj.top.copy()

        print ("_f_ndarray")
        _f_ndarray(f0, xyz)
        print (f0.xyz.shape)

        print ("_f_dontknow")
        f1 = Trajectory()
        f1.top = traj.top.copy()
        _f_dontknow(f1, xyz)
        print (f0.xyz.shape, f1.xyz.shape)
        #assert_almost_equal(f0.xyz, f1.xyz)

if __name__ == "__main__":
    unittest.main()
