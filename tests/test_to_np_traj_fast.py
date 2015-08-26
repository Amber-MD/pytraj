from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        # trp-cage
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        t0 = traj._to_np_traj_fast(0, 10, 1)
        aa_eq(traj.xyz, t0.xyz)
        aa_eq(traj.unitcells, t0.unitcells)

        # tz2
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        t0 = traj._to_np_traj_fast(0, 10, 1)
        aa_eq(traj.xyz, t0.xyz)
        aa_eq(traj.unitcells, t0.unitcells)

if __name__ == "__main__":
    unittest.main()
