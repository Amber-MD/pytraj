from __future__ import print_function
import unittest; import pytraj as pt
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


if __name__ == "__main__":
    unittest.main()
