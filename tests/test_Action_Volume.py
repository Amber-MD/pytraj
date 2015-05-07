import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        dslist = DataSetList()
        act = adict['volume']
        act("", traj, dslist=dslist)
        print (dslist.size)
        print (dslist[0][:])

    def test_0(self):
        from pytraj.common_actions import calc_volume
        traj = mdio.iterload("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")

        # calucate volume for 0-th to 8-th frame, skiip every 2 frames)
        d0 = calc_volume(traj(0, 8, 2), top=traj.top)
        print (d0)

if __name__ == "__main__":
    unittest.main()
