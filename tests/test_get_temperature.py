import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        top = Topology("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        traj = mdio.load("./data/Test_RemdTraj/rem.nc.000", top)
        print (traj.temperature_set)

        farray = traj[:]
        print (farray.temperature_set)

if __name__ == "__main__":
    unittest.main()
