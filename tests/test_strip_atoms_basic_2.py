from __future__ import print_function
import unittest
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        # Aim: simple test for segfault
        # Why?
        # I got segfault when itering `traj` and stripping atoms
        # no error with single frame (f0)
        trajiter = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # 
        f0 = trajiter[0]
        f0.strip_atoms('@CA', trajiter.top, update_top=False)

        print ("start iterating")
        for frame in trajiter:
            frame.strip_atoms('@CA', trajiter.top, update_top=False)

if __name__ == "__main__":
    unittest.main()
