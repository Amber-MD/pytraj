from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        arr0 = np.asarray(f0._buffer2d)
        arr1 = np.asarray(f0)
        aa_eq(arr0.flatten(), arr1.flatten())
        arr0[0, 0] = 100.
        aa_eq(arr0.flatten(), arr1.flatten())
        aa_eq(arr0.flatten(), f0.xyz.flatten())


if __name__ == "__main__":
    unittest.main()
