from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        f0 = traj[0]
        arr0 = np.asarray(f0._buffer2d)
        arr1 = np.asarray(f0)
        aa_eq(arr0.flatten(), arr1.flatten())
        arr0[0, 0] = 100.
        aa_eq(arr0.flatten(), arr1.flatten())
        aa_eq(arr0.flatten(), f0.xyz.flatten())


if __name__ == "__main__":
    unittest.main()
