from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca
from pytraj.testing import Timer


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.datasets import DataSet_Coords_TRJ as DTRJ
        from timeit import timeit
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dtrj = DTRJ()
        dtrj.top = traj.top.copy()
        dtrj.load(traj.filename)

        #print(traj, dtrj)

        #@Timer()
        def test_slice_traj(traj, idx):
            def test():
                traj[idx]

            #print(timeit(test, number=100))

        def test_iter(traj):
            def test():
                for f in traj:
                    pass

            #print(timeit(test, number=100))

        from pytraj.common_actions import calc_radgyr

        def test_action(traj):
            def test():
                calc_radgyr(traj)

            #print(timeit(test, number=100))

            #print("test_slice_traj")
            #print("traj")

        test_slice_traj(traj, 9)
        #print("dtrj")
        test_slice_traj(dtrj, 9)

        #print("test_iter")
        #print("traj")
        test_iter(traj)
        #print("dtrj")
        test_iter(dtrj)

        #print("test_action")
        #print("traj")
        test_action(traj)
        #print("dtrj")
        test_action(dtrj)


if __name__ == "__main__":
    unittest.main()
