from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        from pytraj.utils import convert as cv
        arange = cv.array_to_cpptraj_range(range(7))
        print(arange)

        a0 = pt.multidihedral(traj, resrange='1-7').values
        a1 = pt.multidihedral(traj, resrange=range(7)).values
        aa_eq(a0.flatten(), a1.flatten())


if __name__ == "__main__":
    unittest.main()
