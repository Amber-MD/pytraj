from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    @test_if_having("h5py")
    def test_0(self):
        from pytraj.trajs.TrajectoryHDF5Iterator import TrajectoryHDF5Iterator
        t = TrajectoryHDF5Iterator("./data/ala2.h5")
        print (t.xyz)
        print (t)

        for frame in t:
            print (frame)

if __name__ == "__main__":
    unittest.main()
