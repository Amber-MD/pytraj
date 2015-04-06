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
            f0 = traj[0]
            arr0 = f0.xyz
            print (f0[0])
            arr0[0] = np.array([1, 2, 3])
            print (f0[0])
            print (arr0[0])
            assert_almost_equal(f0[0], arr0[0])
        else:
            print ("need numpy. skip test")

if __name__ == "__main__":
    unittest.main()
