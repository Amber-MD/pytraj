from __future__ import print_function
import unittest
from pytraj.testing import is_linux, test_if_having, aa_eq
from pytraj import io


class Test(unittest.TestCase):

    def test_0(self):
        if is_linux():
            from pytraj.misc import file_type_info
            fname0 = "data/Tc5b.top"
            print(file_type_info(fname0).decode())
        else:
            print("only do this test for Linux")

    @test_if_having("h5py")
    def test_1(self):
        traj = io.load_hdf5("./data/ala2.h5")

        # blindly load file and make guess
        traj2 = io.load("./data/ala2.h5")
        print(traj2)
        aa_eq(traj.xyz, traj2.xyz)

        traj3 = io.load("./data/ala2_h5_faketop.top")
        print(traj3)
        aa_eq(traj.xyz, traj3.xyz)

        traj4 = io.load("./data/ala2.crazy_ext")
        print(traj4)
        aa_eq(traj.xyz, traj4.xyz)

if __name__ == "__main__":
    unittest.main()
