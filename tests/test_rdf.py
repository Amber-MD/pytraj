from __future__ import print_function
import unittest
import numpy as np
import os
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestRDF(unittest.TestCase):
    def test_rdf(self):
        traj = pt.iterload(
            "./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")

        dslist = pt.rdf(traj, "0.5 10.0 :5@CD :WAT@O")

        saved_data_dir = os.path.join(
            cpptraj_test_dir, "Test_Radial", "Radial.agr.save")
        saved_data = np.loadtxt(saved_data_dir, skiprows=8, usecols=(1, ))

        aa_eq(dslist[0], saved_data)
        print(dslist[0])

if __name__ == "__main__":
    unittest.main()
