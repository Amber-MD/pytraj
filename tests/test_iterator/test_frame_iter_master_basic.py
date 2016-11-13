from __future__ import print_function
import unittest
from pytraj import io as mdio


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        indices = traj.top("@CA").indices
        for frame in traj(mask=indices):
            pass

        for frame in traj(mask=[1, 5, 8, 10]):
            pass


if __name__ == "__main__":
    unittest.main()
