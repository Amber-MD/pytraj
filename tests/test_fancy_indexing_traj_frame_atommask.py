from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        aa_eq(traj[0, '@CA'].xyz, traj['@CA'][0].xyz)


if __name__ == "__main__":
    unittest.main()
