# turn off unittest to pass travis
# still need to fix this
import unittest
import numpy as np
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['rmsd']
        dslist = DataSetList()
        act(":3-18@CA reftraj ./data/Tc5b.crd ./data/Tc5b.top", traj, dslist=dslist)
        d0 = dslist[0]
        print (d0[:])
        cpptraj_rmsd = np.loadtxt("./data/rmsd_CA_to_nat_res3_18.dat",
                                  skiprows=1).transpose()[1]
        first10_rmsd = cpptraj_rmsd[:10]
        print(first10_rmsd)

        # make sure to reprofuce cpptraj
        assert_almost_equal(d0[:], first10_rmsd)

if __name__ == "__main__":
    unittest.main()
