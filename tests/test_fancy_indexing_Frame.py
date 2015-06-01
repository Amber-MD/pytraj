from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        xyz = f0.xyz
        atm = traj.top.select("@CA")
        indices = atm.indices
        mask = '@CA'
        f0.set_top(traj.top)
        aa_eq(f0[mask], f0[atm])
        aa_eq(f0[mask], xyz[indices])

        aa_eq(f0['@CA', 0], xyz[indices][0])
        aa_eq(f0['@CA', 0, 0], xyz[indices][0, 0])
        aa_eq(f0[atm, 0], xyz[indices][0])
        aa_eq(f0[atm, 0, 0], xyz[indices][0, 0])
        aa_eq(f0[indices][0, 0], xyz[indices][ 0, 0])

if __name__ == "__main__":
    unittest.main()
