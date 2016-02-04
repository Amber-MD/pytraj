from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir, duplicate_traj


class Test(unittest.TestCase):

    def test_0(self):
        # TODO: 1D or 2D matrix?
        import numpy as np
        from pytraj import matrix as ma
        traj = io.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        saved_data = np.loadtxt("./data/tc5b.matrix_CA.dat")

        arr0 = ma.dist(traj, '@CA', dtype='dataset').to_dict()['Mat_00000']
        dslist1 = ma.dist(traj, '@CA', dtype='dataset')
        arr1 = dslist1.to_dict()['Mat_00000']
        arr2 = ma.dist(traj, '@CA', dtype='dict')['Mat_00000']
        arr3 = ma.dist(traj, '@CA', dtype='dataset').to_ndarray()

        aa_eq(arr0, arr1)
        aa_eq(arr1, arr2.flatten(), decimal=3)
        aa_eq(arr1, arr3.flatten(), decimal=3)
        aa_eq(arr1, saved_data.flatten(), decimal=3)


if __name__ == "__main__":
    unittest.main()
