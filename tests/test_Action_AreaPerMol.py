import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        act = adict['areapermol']
        dslist = DataSetList()
        act("", traj, dslist=dslist)
        print (dslist.size)
        print (dslist[0][:])

if __name__ == "__main__":
    unittest.main()
