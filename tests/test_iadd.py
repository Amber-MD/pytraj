from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.utils import Timer

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        old_size = traj.n_frames
        print (traj)

        traj += traj
        print (traj)
        assert  traj.n_frames == 2 * old_size
        aa_eq(traj[0].xyz, traj[old_size].xyz)

        @Timer()
        def normal(fa):
            fa += fa

        @Timer()
        def make_copy(fa):
            fa += fa.copy()

        fa = traj.copy()
        normal(fa)
        fa = traj.copy()
        make_copy(fa)

        # result: similiar speed

if __name__ == "__main__":
    unittest.main()
