from __future__ import print_function
import unittest
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        print (top[0])
        print (top[-1])
        print (top.atomlist[-1])

if __name__ == "__main__":
    unittest.main()
