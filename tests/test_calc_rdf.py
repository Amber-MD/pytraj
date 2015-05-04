from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    @test_if_path_exists(cpptraj_test_dir)
    def test_0(self):
        import numpy as np
        import os
        traj = mdio.load("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        # calc_rdf with  spacing=0.5, max=10.0, solute mask=":5@CD
        # solvent mask = ":WAT@O"
        dslist = pyca.calc_rdf(traj, "Radial.agr 0.5 10.0 :5@CD :WAT@O")
        print (dslist.size)
        print (dslist[0].tolist())
        # you can find the test here too:
        # https://github.com/mojyt/cpptraj/tree/master/test/Test_Radial/Radial.arg.save
        saved_data_dir = os.path.join(cpptraj_test_dir, "Test_Radial", "Radial.agr.save")
        saved_data = np.loadtxt(saved_data_dir, skiprows=8, usecols=(1,))
        aa_eq(dslist[0].tolist(), saved_data)
        print (dslist.size)

if __name__ == "__main__":
    unittest.main()
