from __future__ import print_function
import unittest
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        # Aim: simple test for segfault
        # bugs : if slicing `traj` (iterator), stripping (_fast_strip_atoms)
        # and slicing again will get segmentation fault
        # *** Error in `python': double free or corruption (!prev):`
        # see (`is slice` in __getitem__ in Trajin.pyx)
        # STATUS: FIXED
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa1 = traj[:]
        fa1._fast_strip_atoms('@CA')
        fa2 = traj[0:10:2]
        fa2._fast_strip_atoms('@CA')
        fa3 = traj[:]
        fa3._fast_strip_atoms('@CA')

        fa3 = traj._fast_slice(slice(2, 10, 1))
        fa3._fast_strip_atoms('@CA')

        fa3 = traj._fast_slice(slice(2, 10, 1))
        fa3._fast_strip_atoms('@CA')

if __name__ == "__main__":
    unittest.main()
