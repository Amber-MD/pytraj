from __future__ import print_function
import pytraj as pt
import unittest; import pytraj as pt
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
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        # use "info("closest") to see its doc (from pytraj import info;
        # info("closest"))
        fa = pyca.closest(traj, ":2,4 center", n_solvents=100)
        pdb_file = cpptraj_test_dir + "/Test_Closest/center.closest.pdb.save"
        saved_frame = mdio.iterload(pdb_file, pdb_file)[0]

        # cpptraj did test for 5-th frame (index starts from 1)
        aa_eq(fa[4].coords, saved_frame.coords, decimal=1)

        fa2, dslist2 = pyca.closest(
            traj, mask=":2,4 center", n_solvents=15,
            restype='all')

if __name__ == "__main__":
    unittest.main()
