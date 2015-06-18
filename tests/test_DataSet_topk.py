from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.testing import aa_eq

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = traj.calc_multidihedral()
        d0 = dslist[0]
        top5 = d0.topk(5)
        arr0 = sorted(d0.values)[-5:-1]
        print (top5, arr0)
        aa_eq(sorted(top5), arr0)

if __name__ == "__main__":
    unittest.main()
