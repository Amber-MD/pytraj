from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload(fn('Test_NAstruct/adh026.3.pdb'))
        d = pt.calc_delta(traj, resrange='1').values

        d1 = pt.dihedral(traj, ":1@C5' :1@C4' :1@C3' :1@O3'")
        d2 = pt._dihedral_res(traj, ("C5'", "C4'", "C3'", "O3'"))

        aa_eq(d, d1)
        aa_eq(d, d2)


if __name__ == "__main__":
    unittest.main()
