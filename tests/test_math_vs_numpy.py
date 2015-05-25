from __future__ import print_function
import unittest
from pytraj.utils import eq, aa_eq, Timer
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import make_fake_traj

class Test(unittest.TestCase):
    def test_0(self):
        @Timer()
        def test_iadd(f0, f1):
            f0.__iadd__(f1)

        @Timer()
        def test_numpy(xyz, xyz0):
            xyz += xyz0

        def test_speed_frame():
            print ("test_speed_frame")
            for n_atoms in (10000, 200000, 400000):
                f00 = make_fake_traj(1, n_atoms)[0]
                xyz = f00.xyz.copy()
                f1 = f00 - 1.
                xyz0 = f1.xyz.copy()
                print ("n_atoms = ", f00.n_atoms)
                print ("test_iadd")
                test_iadd(f00, f1)
                print ("test_numpy")
                test_numpy(xyz, xyz0)
                aa_eq(f00.xyz, xyz)
                print ("")

        def test_speed_traj():
            print ("test_speed_traj")
            for n_atoms in (10000, 200000, 400000):
                f00 = make_fake_traj(10, n_atoms)
                xyz = f00.xyz.copy()
                f1 = f00.copy()
                f1 -= 1.
                xyz0 = f1.xyz.copy()
                print ("n_atoms = ", f00.n_atoms)
                print ("test_iadd")
                test_iadd(f00, f1)
                print ("test_numpy")
                test_numpy(xyz, xyz0)
                aa_eq(f00.xyz, xyz)
                print ("")

        test_speed_traj()


if __name__ == "__main__":
    unittest.main()
