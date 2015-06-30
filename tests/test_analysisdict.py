from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import analdict
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        print(analdict.keys())

    def test_1(self):
        pass
        #traj = mdio.iterload("./data/tz2.nc",  "./data/tz2.parm7")
        #print (traj)
        #analdict['rms2d'](":2@CA :10@CA", traj, traj.top)

if __name__ == "__main__":
    unittest.main()
