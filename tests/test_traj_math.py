from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.utils import Timer

class Test(unittest.TestCase):
    def test_0(self):
        # math with frame_iter
        trajiter = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        saved_xyz = traj.xyz[:].copy()

        traj_view = traj[0:3]
        traj_view += trajiter(stop=3)
        aa_eq(saved_xyz[0:3] * 2, traj[0:3].xyz)

        # reload
        traj = trajiter[:]
        saved_xyz = trajiter.xyz[:].copy()
        traj *= trajiter
        aa_eq(saved_xyz**2, traj.xyz)

    def test_1(self):
        import numpy as np
        #traj = mdio.load("./data/nogit/tip3p/md.trj", "./data/nogit/tip3p/tc5bwat.top")[:5]
        trajiter = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj_saved = trajiter[:]
        traj_saved.join([trajiter[:] for _ in range(200)], copy=False)

        traj = traj_saved.copy()
        print (traj)
        xyz = traj.xyz[:]
        xyz0 = xyz[0].copy()

        # iadd
        traj = traj_saved.copy()
        xyz = traj.xyz[:]
        traj += 1.0
        xyz += 1.0
        aa_eq(traj.xyz, xyz)
        traj += xyz0
        xyz += xyz0
        aa_eq(traj.xyz, xyz)
        
        trajcp = traj.copy()
        xyz_s = traj.xyz.copy()
        trajcp += 2.
        traj += trajcp
        aa_eq(traj.xyz, trajcp.xyz + xyz_s)

        # isub
        traj = traj_saved.copy()
        xyz = traj.xyz[:]
        traj -= 1.0
        xyz -= 1.0
        aa_eq(traj.xyz, xyz)
        fa = traj.copy()
        traj -= fa
        aa_eq(traj.xyz, fa.xyz - fa.xyz)

        # idiv
        traj = traj_saved.copy()
        xyz = traj.xyz[:]
        traj /= 2.0
        xyz /= 2.0
        aa_eq(traj.xyz, xyz)

        xyz_s = traj.xyz.copy()
        traj2 =  traj.copy()
        traj2 /= 0.5
        xyz_2 = traj2.xyz.copy()
        traj /= traj2
        aa_eq(traj.xyz, xyz_s / xyz_2)

        # imul
        traj = traj_saved.copy()
        f0 = traj[0]
        f0 += 1.
        xyz0 = f0.xyz.copy()
        xyz = traj.xyz[:]
        traj *= 2.0
        xyz *= 2.0
        aa_eq(traj.xyz, xyz)

        fa = traj.copy()
        traj *= fa
        aa_eq(traj.xyz, fa.xyz[:]**2)

        @Timer()
        def time_traj(traj, f0=f0):
            traj += 1.
            traj *= 1.
            traj /= 1.
            traj -= 1.
            traj += f0

        @Timer()
        def time_np(xyz, xyz0=xyz0):
            xyz += 1.
            xyz *= 1.
            xyz /= 1.
            xyz -= 1.
            xyz += xyz0

        print ("time_traj")
        time_traj(traj)
        print ("time_np")
        time_np(xyz)

    def test_2(self):
        import numpy as np
        from pytraj.testing import duplicate_traj
        trajiter = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj = duplicate_traj(trajiter, 10000)
        xyz = traj.xyz[:].copy()

        @Timer()
        def test_pytraj_openmp(traj):
            traj += traj

        @Timer()
        def test_numpy(xyz):
            xyz += xyz

        print ("test_pytraj_openmp")
        test_pytraj_openmp(traj)
        print ("test_numpy")
        test_numpy(xyz)

if __name__ == "__main__":
    unittest.main()
