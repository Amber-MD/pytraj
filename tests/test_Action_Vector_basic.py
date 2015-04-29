from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import test_if_having

class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        from pytraj import calculate
        traj = mdio.load("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        dslist = calculate("vector", traj, "@CA @CB mass")
        print ('dslist.size = ', dslist.size)
        print (dslist.get_dtypes())
        for d0 in dslist:
            print (d0)
            print (d0.to_ndarray().shape)

    @test_if_having("numpy")
    def test_1(self):
        from pytraj.common_actions import calc_vector
        traj = mdio.load("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        d0 = calc_vector(traj, "@CA @CB mass")
        print (d0.to_ndarray())
        print (d0.tolist())

if __name__ == "__main__":
    unittest.main()
