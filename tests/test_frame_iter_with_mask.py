from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        traj0_CA = traj[:]
        traj0_CA.strip_atoms("!@CA")

        # test TrajReadOnly
        for idx, f0 in enumerate(traj(mask='@CA')):
            f1 = traj0_CA[idx]
            print (idx, f0, f1)
            assert_almost_equal(f0.coords, f1.coords)
            print (f0.rmsd(f1))

        # test FrameArray
        for idx, f0 in enumerate(farray(mask='@CA')):
            f1 = traj0_CA[idx]
            print (idx, f0, f1)
            assert_almost_equal(f0.coords, f1.coords)
            print (f0.rmsd(f1))

if __name__ == "__main__":
    unittest.main()
