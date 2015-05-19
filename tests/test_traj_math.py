from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.utils import Timer

class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        traj = mdio.load("./data/nogit/tip3p/md.trj", "./data/nogit/tip3p/tc5bwat.top")[:5]
        print (traj)
        xyz = traj.xyz[:]

        # iadd
        traj += 1.0
        xyz += 1.0
        aa_eq(traj.xyz, xyz)

        # isub
        traj -= 1.0
        xyz -= 1.0
        aa_eq(traj.xyz, xyz)

        # idiv
        traj /= 2.0
        xyz /= 2.0
        aa_eq(traj.xyz, xyz)

        # imul
        traj *= 2.0
        xyz *= 2.0
        aa_eq(traj.xyz, xyz)

        @Timer()
        def time_traj(traj):
            traj += 1.
            traj *= 1.
            traj /= 1.
            traj -= 1.

        @Timer()
        def time_np(xyz):
            xyz += 1.
            xyz *= 1.
            xyz /= 1.
            xyz -= 1.

        print ("time_traj")
        time_traj(traj)
        print ("time_np")
        time_np(xyz)

if __name__ == "__main__":
    unittest.main()
