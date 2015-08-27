from __future__ import print_function
import unittest
import pytraj as pt
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


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
        for f0, f1 in zip(traj0(autoimage=True,
                                rmsfit=(3, '@CA,C,N')),
                          traj1(autoimage=True,
                                rmsfit=(3, '@CA,C,N'))):
            assert f0.same_coords_as(f1)
            assert f0.rmsd_nofit(f1) == 0.


if __name__ == "__main__":
    unittest.main()
