import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        f0 = traj[0]
        f0saved = f0.copy()
        f1 = traj[1]

        rmsd_0 = f0.rmsd(f1)
        rmsd_0_nofit = f0.rmsd_nofit(f1)
        assert rmsd_0 != rmsd_0_nofit

        # do fitting
        f1.fit_to(f0)
        rmsd_1 = f1.rmsd(f0)
        rmsd_1_nofit = f1.rmsd_nofit(f0)

        # make sure that rmsd_nofit after do fitting is equal to rmsd (with fit)
        print (rmsd_1, rmsd_1_nofit)
        assert rmsd_1 - rmsd_1_nofit < 1E-3

        print ("test farray.fit_to")
        farray.fit_to(f0)
        print (farray[1].rmsd_nofit(f0))
        print (farray[1].rmsd(f0))
        assert rmsd_1 - farray[1].rmsd_nofit(f0) < 1E-3
        print (farray.fit_to)

if __name__ == "__main__":
    unittest.main()
