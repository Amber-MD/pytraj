from __future__ import print_function
import unittest
from pytraj import *
from pytraj import io as mdio


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        f0 = Frame()
        f0.append_xyz(traj[0]._buffer2d)


if __name__ == "__main__":
    unittest.main()
