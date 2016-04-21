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
        from pytraj import dihedral_analysis as da
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        # resrange 7, phi psi
        saved_file = './data/test_multidihedral.dat'
        saved_data_phi = np.loadtxt(saved_file, skiprows=1).transpose()[1]

        arr0 = da.calc_phi(traj, resrange='7').to_ndarray()
        dslist1 = da.calc_phi(traj, resrange='7')
        arr1 = dslist1.to_ndarray()
        arr2 = da.calc_phi(traj, resrange='7', dtype='ndarray')
        arr3 = da.calc_phi(traj, resrange='7').to_dict()

        # assert to cpptraj output
        aa_eq(arr0, saved_data_phi)
        aa_eq(arr1, saved_data_phi)
        aa_eq(arr2, saved_data_phi)

        for key in dslist1.keys():
            aa_eq(dslist1[key].values, arr3[key])


if __name__ == "__main__":
    unittest.main()
