import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        f0 = traj[0]
        act = adict['autoimage']
        f0cp = f0.copy()
        assert f0.same_coords_as(f0cp) == True
        act("", f0, traj.top)
        assert f0.same_coords_as(f0cp) == False

    def test_1(self):
        traj = mdio.iterload("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        f0 = traj[0]
        f0cp = f0.copy()
        print (f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == True
        adict['autoimage']("", f0, traj.top)
        print (f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == False

        fsaved = mdio.iterload("./data/tz2.truncoct.autoiamge.save.r",
                           "./data/tz2.truncoct.parm7")[0]
        assert_almost_equal(fsaved.coords, f0.coords)

    def test_2(self):
        from pytraj.common_actions import do_autoimage
        # test do_autoimage
        traj = mdio.iterload("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        f0 = traj[0]
        f0cp = f0.copy()
        print (f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == True
        do_autoimage(traj=f0, top=traj.top)
        print (f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == False

        fsaved = mdio.iterload("./data/tz2.truncoct.autoiamge.save.r",
                           "./data/tz2.truncoct.parm7")[0]
        assert_almost_equal(fsaved.coords, f0.coords)

if __name__ == "__main__":
    unittest.main()
