from __future__ import print_function
import unittest
from pytraj import io as mdio


class Test(unittest.TestCase):
    def test_0(self):
        # Aim: simple test for segfault
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        fa.strip_atoms("@CA")


if __name__ == "__main__":
    unittest.main()
