from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.six_2 import zip
from timeit import timeit

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/nogit/remd/remd.x.000", "./data/nogit/remd/myparm.top")[:200]

        s = slice(3, 100, 4)
        t0 = traj[s]
        t1 = traj._fast_slice(s)
        print (traj[3, 0], t0[0, 0], t1[0, 0])
        aa_eq(t0.xyz, t1.xyz)
        assert t0.n_frames == t1.n_frames

        # test timer
        def normal_slice():
            traj[s]

        def fast_slice():
            traj._fast_slice(s)

        print (timeit(normal_slice, number=100))
        print (timeit(fast_slice, number=100))

if __name__ == "__main__":
    unittest.main()
