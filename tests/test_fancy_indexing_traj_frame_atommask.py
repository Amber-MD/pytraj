from __future__ import print_function
import pytraj as pt
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0, '@CA']
        f1 = traj['@CA'][0]
        assert pt.tools.rmsd(f0.xyz, f1.xyz) == 0.0


if __name__ == "__main__":
    unittest.main()
