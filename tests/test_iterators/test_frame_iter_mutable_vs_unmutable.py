from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.compat import zip
        from pytraj import Trajectory, TrajectoryIterator
        # TrajectoryIterator object
        traj0 = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        assert isinstance(traj0, TrajectoryIterator)

        # Trajectory object
        traj1 = traj0[:]
        assert isinstance(traj1, Trajectory)

        # check same coords
        ref = traj0[0]

        for f0, f1 in zip(
                traj0(autoimage=True,
                      rmsfit=(ref, '@CA,C,N')),
                traj1(autoimage=True,
                      rmsfit=(ref, '@CA,C,N'))):
            aa_eq(f0.xyz, f1.xyz)
            assert f0.rmsd_nofit(f1) == 0.


if __name__ == "__main__":
    unittest.main()
