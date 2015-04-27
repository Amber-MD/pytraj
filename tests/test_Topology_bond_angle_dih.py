from __future__ import print_function
import unittest
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top

        # numbers wer taken from top.summary() (C++ stdout)
        assert list(top.bonds).__len__() == 310 
        assert list(top.angles).__len__() == 565 
        assert list(top.dihedrals).__len__() == 1351

if __name__ == "__main__":
    unittest.main()
