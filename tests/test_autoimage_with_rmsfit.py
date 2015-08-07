from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq, eq_coords
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca
from pytraj.misc import rmsd as rmsd_2darray


class Test(unittest.TestCase):
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
            print(rmsd_2darray(f1.xyz, f2.xyz))

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
