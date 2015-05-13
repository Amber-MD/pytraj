from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from timeit import timeit

class Test(unittest.TestCase):
    def test_0(self):
        # Trajectory
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        for i in range(5):
            traj += traj.copy()
        traj2 = traj.copy()
        print (traj, traj2)
        xyz = traj.xyz * 2

        def normal():
            traj.update_xyz(xyz)

        def use_memcpy():
            traj2._update_xyz_memcpy(xyz)
        
        print (timeit(normal, number=1000))
        print (timeit(use_memcpy, number=1000))

        aa_eq(traj.xyz, xyz)
        aa_eq(traj2.xyz, xyz)

    def test_1(self):
        # Frame
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0, f1 = traj[:2]
        xyz = f0.xyz + 1.0
        f0._fast_copy_from_xyz(xyz)
        aa_eq(f0.xyz, xyz)
        f1._fast_copy_from_frame(f0)
        aa_eq(f0.xyz, f1.xyz)

if __name__ == "__main__":
    unittest.main()
