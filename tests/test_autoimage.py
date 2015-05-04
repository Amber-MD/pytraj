import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
import pytraj.io as io
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import test_if_having

class Test(unittest.TestCase):
    def test_1(self):
        traj = mdio.load("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        f0 = traj[0]
        f0cp = f0.copy()
        print (f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == True
        adict['autoimage']("", f0, traj.top)
        print (f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == False

        fsaved = mdio.load("./data/tz2.truncoct.autoiamge.save.r",
                           "./data/tz2.truncoct.parm7")[0]
        assert_almost_equal(fsaved.coords, f0.coords)

    def test_2(self):
        from pytraj.common_actions import do_autoimage
        # test do_autoimage
        traj = mdio.load("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        f0 = traj[0]
        f0cp = f0.copy()
        print (f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == True
        do_autoimage(traj=f0, top=traj.top)
        print (f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == False

        fsaved = mdio.load("./data/tz2.truncoct.autoiamge.save.r",
                           "./data/tz2.truncoct.parm7")[0]
        assert_almost_equal(fsaved.coords, f0.coords)

    @test_if_having("mdtraj")
    def test_3(self):
        import mdtraj as md
        import numpy as np
        # load from mdtraj
        # test do_autoimage
        traj = mdio.load("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        m_traj = md.load_netcdf("./data/tz2.truncoct.nc", top="./data/tz2.truncoct.parm7")
        # create mutable FrameArray from TrajReadOnly
        fa = traj[:]
        # do autoiamge
        fa.autoimage()

        # try loading from `mdtraj` object
        fake_fa = io.load_mdtraj(m_traj)
        # try autoiamging
        fake_fa.autoimage()
        assert_almost_equal(fa.xyz.flatten(), fake_fa.xyz.flatten(), decimal=3)

        for frame in fake_fa:
            assert frame.has_box() == True
        assert fake_fa.top.has_box() == True

if __name__ == "__main__":
    unittest.main()
