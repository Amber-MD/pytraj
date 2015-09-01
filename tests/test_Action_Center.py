from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        fa = traj[:]

        saved_traj = pt.load("./data/tz2.center_mass.nc", traj.top)
        fa.center(":1 mass")
        aa_eq(fa[:2].xyz, saved_traj.xyz, decimal=5)


if __name__ == "__main__":
    unittest.main()
