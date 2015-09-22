from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.testing import make_fake_traj, Timern

class Test(unittest.TestCase):
    def test_0(self):
        def test_copy(n_atoms):
            # test copy
            # Conclusion: similiar speed between numpy and _fast_copy_from_xyz
            print ("-------------------------------------")
            print ('n_atoms', n_atoms)
            frame = make_fake_traj(1, n_atoms)[0]
            xyz = frame.xyz.copy()
            xyz0 = frame.xyz.copy()

            @Timer()
            def test_pytraj_fast_copy(frame, xyz0):
                frame._fast_copy_from_xyz(xyz0)

            @Timer()
            def test_pytraj_xyz(frame, xyz0):
                frame.xyz[:] = xyz0

            @Timer()
            def test_numpy(xyz, xyz0):
                xyz[:] = xyz0

            print ("test_pytraj_fast_copy")
            test_pytraj_fast_copy(frame, xyz0)
            print ("")
            print ("test_pytraj_xyz")
            test_pytraj_xyz(frame, xyz0)
            print ("")
            print ("test_numpy")
            test_numpy(xyz, xyz0)
            print ("")
            aa_eq(frame.xyz, xyz)

        for n_atoms in (100, 1000, 10000, 100000, 500000):
            test_copy(n_atoms)

if __name__ == "__main__":
    unittest.main()
