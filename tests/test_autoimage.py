import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj.compat import zip
from pytraj.testing import aa_eq


class TestRegular(unittest.TestCase):
    def test_1(self):
        traj = pt.iterload(
            "./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        f0 = traj[0]
        f0cp = f0.copy()
        #print(f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == True
        adict['autoimage']("", f0, traj.top)
        #print(f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == False

        fsaved = pt.iterload("./data/tz2.truncoct.autoiamge.save.r",
                             "./data/tz2.truncoct.parm7")[0]
        aa_eq(fsaved.coords, f0.coords, decimal=3)

    def test_2(self):
        from pytraj.common_actions import do_autoimage
        # test do_autoimage
        traj = pt.iterload(
            "./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        f0 = traj[0]
        f0cp = f0.copy()
        #print(f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == True
        do_autoimage(traj=f0, top=traj.top)
        #print(f0.same_coords_as(f0cp))
        assert f0.same_coords_as(f0cp) == False

        fsaved = pt.iterload("./data/tz2.truncoct.autoiamge.save.r",
                             "./data/tz2.truncoct.parm7")[0]
        aa_eq(fsaved.coords, f0.coords, decimal=3)

    def test_4(self):
        # combined with get_coordinates
        from pytraj.misc import rmsd_1darray
        traj0 = pt.iterload(
            "./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")

        # test autoimage
        traj1 = traj0[:]
        xyz0 = pt.get_coordinates(traj0(autoimage=True))
        traj1.autoimage()
        xyz1 = traj1.xyz
        # OK
        assert rmsd_1darray(xyz0.flatten(), xyz1.flatten()) == 0.

        # rmsfit
        # reset traj1
        traj1 = traj0[:]
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
        traj1 = traj0[:]
        # get new trajectory from traj0
        traj2 = pt._load_from_frame_iter(
            traj0(autoimage=True,
                  rmsfit=(0, '@CA,C,N')))
        traj1.autoimage()
        traj1.rmsfit(ref=0, mask='@CA,C,N')

        # take '@CA,C,N' coords
        xyz2 = traj2['@CA,C,N'].xyz
        xyz1 = traj1['@CA,C,N'].xyz
        # PASSED
        assert rmsd_1darray(xyz1.flatten(), xyz2.flatten()) < 1E-10


class TestWithRmsfit(unittest.TestCase):
    def test_0(self):
        # TrajectoryIterrator
        # status: failed
        from pytraj.compat import zip
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        pt.write_traj("./output/tz2.autoimage_with_rmsfit.nc",
                      traj(autoimage=True,
                           rmsfit=(0, '@CA,C,N')),
                      overwrite=True)

        saved_traj = pt.load('data/tz2.autoimage_with_rmsfit.nc', traj.top)
        p_traj = pt.load('./output/tz2.autoimage_with_rmsfit.nc', traj.top)

        aa_eq(saved_traj.xyz, p_traj.xyz)
        for f1, f2 in zip(p_traj, saved_traj):
            pass

    def test_1(self):
        # Trajectory (mutable)
        # status: OK
        from pytraj.compat import zip
        # note: use `load` instead of `iterload`
        traj = pt.load("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        traj.autoimage()
        traj.rmsfit(mask='@CA,C,N')
        saved_traj = pt.load('data/tz2.autoimage_with_rmsfit.nc', traj.top)

        # PASSED
        aa_eq(saved_traj.xyz, traj.xyz)


if __name__ == "__main__":
    unittest.main()
