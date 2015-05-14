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
from pytraj.externals.six.moves import range
from timeit import timeit

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")[:]
        print (traj.n_frames)

        s = slice(3, traj.n_frames, 4)
        t0 = traj[s]
        t1 = traj._fast_slice(s)
        print (traj[3, 0], t0[0, 0], t1[0, 0])
        aa_eq(t0.xyz, t1.xyz)
        assert t0.n_frames == t1.n_frames

        # test timer
        def normal_slice():
            traj[s]

        sr = range(3, traj.n_frames, 4)
        def test_range():
            traj[sr]

        def fast_slice():
            traj._fast_slice(s)

        print ("normal_slice")
        print (timeit(normal_slice, number=100))
        print ("test_range")
        print (timeit(test_range, number=100))
        print ("fast_slice")
        print (timeit(fast_slice, number=100))

        #
        aa_eq(traj[sr].xyz, traj[s].xyz)
        aa_eq(traj._fast_slice(s).xyz, traj[s].xyz)

if __name__ == "__main__":
    unittest.main()
