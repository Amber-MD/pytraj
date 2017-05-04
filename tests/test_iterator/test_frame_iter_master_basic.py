from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        indices = traj.top("@CA").indices
        for frame in traj(mask=indices):
            pass

        for frame in traj(mask=[1, 5, 8, 10]):
            pass


if __name__ == "__main__":
    unittest.main()
