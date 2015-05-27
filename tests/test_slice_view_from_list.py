from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        traj = io.load_sample_data("tz2")[:]
        mylist = [1, 5, 8]

        # slicing
        fa0 = traj[mylist]
        fa1 = traj[mylist]
        aa_eq(fa0.xyz, fa1.xyz)
        aa_eq(fa0.xyz, traj[mylist].xyz)

        # update fa0 and make sure fa1, traj are updated too
        fa0[0, 0] = [100., 101., 102.]
        assert fa1[0, 0, 0] == fa0[0, 0, 0] == 100.
        assert traj[1, 0, 0] == 100.
        print (fa0[0 , 0], fa1[0 , 0], traj[1 , 0])

if __name__ == "__main__":
    unittest.main()
