import os
import unittest
import pytraj as pt
from utils import fn

from pytraj import adict
from pytraj.testing import aa_eq, cpptraj_test_dir
from pytraj.externals.six import zip
from pytraj.utils.tools import rmsd_1darray


class TestRegular(unittest.TestCase):
    def test_1(self):
        traj = pt.iterload(fn('tz2.truncoct.nc'), fn('tz2.truncoct.parm7'))
        f0 = traj[0]
        f0.copy()
        adict['autoimage']("", f0, traj.top)

        fsaved = pt.iterload(
            fn('tz2.truncoct.autoiamge.save.r'), fn('tz2.truncoct.parm7'))[0]
        aa_eq(fsaved.xyz, f0.xyz, decimal=3)

    def test_2(self):
        from pytraj.all_actions import do_autoimage
        # test do_autoimage
        traj = pt.iterload(fn('tz2.truncoct.nc'), fn('tz2.truncoct.parm7'))
        f0 = traj[0]
        f0.copy()
        do_autoimage(traj=f0, top=traj.top)

        fsaved = pt.iterload(
            fn('tz2.truncoct.autoiamge.save.r'), fn('tz2.truncoct.parm7'))[0]
        aa_eq(fsaved.xyz, f0.xyz, decimal=3)

    def test_4(self):
        # combined with get_coordinates
        traj0 = pt.iterload(fn('tz2.truncoct.nc'), fn('tz2.truncoct.parm7'))

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
        traj2 = pt.load_from_frame_iter(traj0(rmsfit=(0, '@CA,C,N')))
        traj1.rmsfit(ref=0, mask='@CA,C,N')

        # take '@CA,C,N' xyz
        xyz2 = traj2['@CA,C,N'].xyz
        xyz1 = traj1['@CA,C,N'].xyz
        # OK
        assert rmsd_1darray(xyz1.flatten(), xyz2.flatten()) < 1E-10

        # combine autoimage with rmsfit
        # reset traj1
        traj1 = traj0[:]
        # get new trajectory from traj0
        traj2 = pt.load_from_frame_iter(
            traj0(autoimage=True, rmsfit=(0, '@CA,C,N')))
        traj1.autoimage()
        traj1.rmsfit(ref=0, mask='@CA,C,N')

        # take '@CA,C,N' xyz
        xyz2 = traj2['@CA,C,N'].xyz
        xyz1 = traj1['@CA,C,N'].xyz
        # PASSED
        assert rmsd_1darray(xyz1.flatten(), xyz2.flatten()) < 1E-10


class TestWithRmsfit(unittest.TestCase):
    def test_on_disk_trajectory(self):
        # TrajectoryIterrator
        output = "ok_to_delete.nc"
        traj = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
        pt.write_traj(
            output,
            traj(autoimage=True, rmsfit=(0, '@CA,C,N')),
            overwrite=True)

        saved_traj = pt.load(fn('tz2.autoimage_with_rmsfit.nc'), traj.top)
        p_traj = pt.load(output, traj.top)

        aa_eq(saved_traj.xyz, p_traj.xyz)
        # ensure iterable without memory freeing
        for f1, f2 in zip(p_traj, saved_traj):
            aa_eq(f1.xyz, f2.xyz)

    def test_in_memory_trajectory(self):
        traj = pt.load(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
        traj.autoimage()
        traj.rmsfit(mask='@CA,C,N')
        saved_traj = pt.load(fn('tz2.autoimage_with_rmsfit.nc'), traj.top)
        aa_eq(saved_traj.xyz, traj.xyz)


def test_autoimage_for_tightly_packed_systems():
    trajin_fn = os.path.join(cpptraj_test_dir, 'Test_AutoImage', 'G3_3A.rst7')
    prmtop_fn = os.path.join(cpptraj_test_dir, 'Test_AutoImage',
                             'nowat.G3_3A.parm7')
    saved_rst7 = os.path.join(cpptraj_test_dir, 'Test_AutoImage',
                              'image.G3_3A.rst7.save')
    saved_traj = pt.iterload(saved_rst7, prmtop_fn)

    traj = pt.load(trajin_fn, prmtop_fn)
    pt.autoimage(traj, 'anchor :96 origin')
    aa_eq(traj.xyz, saved_traj.xyz)


if __name__ == "__main__":
    unittest.main()
