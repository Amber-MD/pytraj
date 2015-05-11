import unittest
import numpy as np
from pytraj.base import *
from pytraj import io as io
from pytraj.utils import aa_eq
from pytraj.decorators import test_if_having

class Test(unittest.TestCase):
    @test_if_having("h5py")
    def test_1(self):
        traj = io.load_hdf5("./data/ala2.h5")

        # blindly load file and make guess
        traj2 = io.load("./data/ala2.h5")
        print (traj2)
        aa_eq(traj.xyz, traj2.xyz)

        traj3 = io.load("./data/ala2_h5_faketop.top")
        print (traj3)
        aa_eq(traj.xyz, traj3.xyz)

        traj4 = io.load("./data/ala2.crazy_ext")
        print (traj4)
        aa_eq(traj.xyz, traj4.xyz)

if __name__ == "__main__":
    unittest.main()
