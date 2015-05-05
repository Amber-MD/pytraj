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
        # test dir: 
        # 1. $AMBERHOME/AmberTools/test/cpptraj/Test_Closest/
        # 2. https://github.com/mojyt/cpptraj/tree/master/test/Test_Closest

        import numpy as np
        traj = mdio.load("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        print (traj.top)
        # use "info("closest") to see its doc (from pytraj import info; info("closest"))
        fa = pyca.closest(traj, "10 :2,4 center")
        pdb_file = cpptraj_test_dir + "/Test_Closest/center.closest.pdb.save"
        saved_frame = mdio.load(pdb_file, pdb_file)[0]

        # cpptraj did test for 5-th frame (index starts from 1)
        print (fa[4].rmsd(saved_frame))
        aa_eq(fa[4].coords, saved_frame.coords, decimal=1)

if __name__ == "__main__":
    unittest.main()
