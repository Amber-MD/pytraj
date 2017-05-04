from __future__ import print_function
import unittest
from pytraj import *
import pytraj as pt
from utils import fn


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        f0 = Frame()
        f0.append_xyz(traj[0]._buffer2d)


if __name__ == "__main__":
    unittest.main()
