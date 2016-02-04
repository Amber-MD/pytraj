from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir, duplicate_traj


class Test(unittest.TestCase):

    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        # 3 atoms x 2
        arr = np.array([[12, 13, 15], [16, 17, 18]])

        f0saved = traj[0].copy()
        f0 = traj[0].copy()
        f0.swap_atoms(arr)
        aa_eq(f0saved[12], f0[16])
        aa_eq(f0saved[13], f0[17])
        aa_eq(f0saved[15], f0[18])


if __name__ == "__main__":
    unittest.main()
