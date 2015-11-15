from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestCenter(unittest.TestCase):

    def test_center(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        saved_traj = pt.load("./data/tz2.center_mass.nc", traj.top)

        fa = traj[:]
        pt.center(fa, mask=':1', mass=True)
        aa_eq(fa.xyz, fa.xyz, decimal=5)

        # raise if center not in 'origin', 'box'
        self.assertRaises(ValueError, lambda: pt.center(fa, center='oh'))


if __name__ == "__main__":
    unittest.main()
