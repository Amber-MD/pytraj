from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj[0, '@CA']

if __name__ == "__main__":
    unittest.main()
