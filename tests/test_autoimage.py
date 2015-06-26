import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
import pytraj.io as io
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import test_if_having, no_test
from pytraj.compat import izip as zip

class Test(unittest.TestCase):
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

    @no_test
    @test_if_having("mdtraj")
    def test_3(self):
        import mdtraj as md
        import numpy as np
        # load from mdtraj
        # test do_autoimage
        traj = mdio.iterload("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        fa = traj[:]
        for f0, f1 in zip(traj, fa):
            assert f0.box.type == traj.top.box.type == f1.box.type  == 'truncoct'
        m_traj = md.load_netcdf("./data/tz2.truncoct.nc", top="./data/tz2.truncoct.parm7")
        # create mutable Trajectory from TrajectoryIterator
        fake_fa = io.load_mdtraj(m_traj)
        for frame in fake_fa:
            assert frame.box.type == 'truncoct'
        assert_almost_equal(fa.xyz.flatten(), traj.xyz.flatten())
        assert_almost_equal(fa.xyz.flatten(), fake_fa.xyz.flatten())
        # do autoiamge
        fa.autoimage()
        print (fa[0, 0])
        # try loading from `mdtraj` object
        fake_fa.autoimage()
        print (fake_fa[0, 0])
        assert_almost_equal(fa.xyz.flatten(), fake_fa.xyz.flatten(), decimal=1)

        for frame in fake_fa:
            assert frame.has_box() == True
        assert fake_fa.top.has_box() == True

    def test_4(self):
        # combined with get_coordinates
        from pytraj.misc import rmsd_1darray
        traj0 = mdio.iterload("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")

        # test autoimage
        traj1 = traj0.to_mutable_trajectory()
        xyz0 = pt.get_coordinates(traj0(autoimage=True))
        traj1.autoimage()
        xyz1 = traj1.xyz
        # OK
        assert rmsd_1darray(xyz0.flatten(), xyz1.flatten()) == 0.

        # rmsfit
        # reset traj1
        traj1 = traj0.to_mutable_trajectory()
        # get new trajectory from traj0
        traj2 = pt._load_from_frame_iter(traj0(rmsfit=(0, '@CA,C,N')))
        traj1.rmsfit(ref=0, mask='@CA,C,N')

        # take '@CA,C,N' coords
        xyz2 = traj2['@CA,C,N'].xyz
        xyz1 = traj1['@CA,C,N'].xyz
        # OK
        assert rmsd_1darray(xyz1.flatten(), xyz2.flatten()) < 1E-10

        # combine autoimage with rmsfit
        # reset traj1
        traj1 = traj0.to_mutable_trajectory()
        # get new trajectory from traj0
        traj2 = pt._load_from_frame_iter(traj0(autoimage=True, rmsfit=(0, '@CA,C,N')))
        traj1.autoimage()
        traj1.rmsfit(ref=0, mask='@CA,C,N')

        # take '@CA,C,N' coords
        xyz2 = traj2['@CA,C,N'].xyz
        xyz1 = traj1['@CA,C,N'].xyz
        # Failed
        print(rmsd_1darray(xyz1.flatten(), xyz2.flatten()))
        assert rmsd_1darray(xyz1.flatten(), xyz2.flatten()) < 1E-10


if __name__ == "__main__":
    unittest.main()
