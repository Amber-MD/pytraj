from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        # TODO: 1D or 2D matrix?
        import numpy as np
        from pytraj import matrix_analysis as ma
        traj = io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        saved_data = np.loadtxt("./data/tc5b.matrix_CA.dat")

        arr0 = ma.distance_matrix(traj, '@CA').to_dict()['Mat_00000']
        dslist1 = ma.distance_matrix(traj, '@CA')
        arr1 = dslist1.to_dict()['Mat_00000']
        arr2 = ma.distance_matrix(traj, '@CA', dtype='dict')['Mat_00000']
        arr3 = ma.distance_matrix(traj, '@CA').to_ndarray()

        aa_eq(arr0, arr1)
        aa_eq(arr1, arr2.flatten())
        aa_eq(arr1, arr3.flatten())
        aa_eq(arr1, saved_data.flatten())

if __name__ == "__main__":
    unittest.main()
