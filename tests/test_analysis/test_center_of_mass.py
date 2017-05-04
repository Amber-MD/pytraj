from __future__ import print_function
import pytraj as pt
from utils import fn
import unittest
from pytraj.testing import aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np

        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        traj2 = traj[:]
        d1 = pt.calc_center_of_mass(traj, dtype='dataset')
        d2 = pt.calc_center_of_mass(traj2, dtype='dataset')

        for frame in traj:
            pass

        saved_d0 = np.loadtxt(fn('vec.out'), skiprows=1, usecols=(1, 2, 3))

        aa_eq(d1.to_ndarray().flatten(), saved_d0.flatten())
        aa_eq(d2.to_ndarray().flatten(), saved_d0.flatten())

        aa_eq(
            pt.center_of_geometry(traj, dtype='ndarray'),
            pt.center_of_geometry(traj2, dtype='ndarray'))


if __name__ == "__main__":
    unittest.main()
