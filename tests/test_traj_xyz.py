from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import has_
from pytraj.utils.check_and_assert import assert_almost_equal, eq

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        if has_("numpy"):
            import numpy as np
            from numpy.testing import assert_almost_equal
            arr0 = traj.xyz
            print (arr0.shape)
            assert_almost_equal(arr0, traj[:, :, :])

            # create FrameArray
            farray = traj[:]
            assert_almost_equal(arr0, farray[:, :, :])
        else:
            print ("need numpy. skip test")

if __name__ == "__main__":
    unittest.main()
