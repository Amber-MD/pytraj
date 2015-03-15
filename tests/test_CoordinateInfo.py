import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        cinfo = traj.top.get_parm_coord_info()
        print (dir(cinfo))
        print (cinfo.has_box())
        print (cinfo.has_time())
        print (cinfo.has_vel())
        print (cinfo.traj_box())

if __name__ == "__main__":
    unittest.main()
