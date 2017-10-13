from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        from pytraj.utils import convert as cv
        cv.array_to_cpptraj_range(range(7))

        a0 = pt.multidihedral(traj, resrange='1-7').values
        a1 = pt.multidihedral(traj, resrange=range(7)).values
        aa_eq(a0.flatten(), a1.flatten())


if __name__ == "__main__":
    unittest.main()
