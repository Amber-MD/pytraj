from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.compat import zip
from pytraj.externals.six.moves import range
from timeit import timeit


class Test(unittest.TestCase):
    #@test_if_path_exists("./data/nogit/tip3p/")

    def test_0(self):
        try:
            tip3pdir = "./data/nogit/tip3p/"
            traj = mdio.iterload(
                tip3pdir + "/md.trj", tip3pdir + "/tc5bwat.top")[:500]
            start, stop, step = 0, 277, 3
        except:
            # load small file
            traj = mdio.load("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
            for _ in range(5):
                traj.join((traj.copy(), traj.copy()))
            start, stop, step = 0, 1000, 4
        #print(traj.n_frames)

        s = slice(start, stop, step)
        #print(s)
        t0 = traj[s]
        t1 = traj._fast_slice(s)
        #print(traj[3, 0], t0[0, 0], t1[0, 0])
        aa_eq(t0.xyz, t1.xyz)
        assert t0.n_frames == t1.n_frames

        # test timer
        def normal_slice():
            traj[s]

        sr = range(start, stop, step)

        def test_range():
            traj[sr]

        import numpy as np
        indices = np.random.randint(0, stop, len(sr))
        from array import array
        indices_p = array('i', indices)
        #print(indices_p)

        def test_numpy_array_slicing():
            traj[indices]

        def test_pyarray_slicing():
            traj[indices_p]

        assert traj[indices].n_frames == traj[sr].n_frames

        def fast_slice():
            traj._fast_slice(s)

        def fast_slice2():
            traj._fast_slice((start, stop, step))

        #print("normal_slice")
        #print(timeit(normal_slice, number=100))
        #print("test_range")
        #print(timeit(test_range, number=100))
        #print("test_numpy_array_slicing")
        #print(timeit(test_numpy_array_slicing, number=100))
        #print("test_pyarray_slicing")
        #print(timeit(test_pyarray_slicing, number=100))
        #print("fast_slice")
        #print(timeit(fast_slice, number=100))
        #print("fast_slice2")
        #print(timeit(fast_slice2, number=100))

        #
        aa_eq(traj[sr].xyz, traj[s].xyz)
        aa_eq(traj._fast_slice(s).xyz, traj[s].xyz)
        aa_eq(traj._fast_slice((start, stop, step)).xyz, traj[s].xyz)
        aa_eq(traj[indices].xyz, traj[indices_p].xyz)


if __name__ == "__main__":
    unittest.main()
