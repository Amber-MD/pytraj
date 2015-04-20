from __future__ import print_function
import unittest
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = traj.search_hbonds()
        print (dslist.groupby("SER").to_dict())

if __name__ == "__main__":
    unittest.main()
