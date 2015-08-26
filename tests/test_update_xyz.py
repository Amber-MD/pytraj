from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from timeit import timeit


class Test(unittest.TestCase):
    def test_0(self):
        # Trajectory
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        for i in range(5):
            traj.join(traj.copy())
        print(traj)
        traj2 = traj.copy()
        traj3 = traj.copy()
        print(traj, traj2)
        xyz = traj.xyz * 2

        def normal():
            traj.update_xyz(xyz)

        print(timeit(normal, number=1000))

        aa_eq(traj.xyz, xyz)

    def test_1(self):
        print("Frame")
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0, f1 = traj[:2]
        xyz = f0.xyz + 1.0
        f0._fast_copy_from_xyz(xyz)
        aa_eq(f0.xyz, xyz)
        f1._fast_copy_from_frame(f0)
        aa_eq(f0.xyz, f1.xyz)

        frame = f0.copy()
        xyz = frame.xyz.copy()
        z = xyz.copy()

        def normal():
            frame.buffer2d[:] = xyz

        frame = f0.copy()
        xyz = frame.xyz.copy()
        z = xyz.copy()

        def use_memcpy():
            frame._fast_copy_from_xyz(xyz)

        frame = f0.copy()
        xyz = frame.xyz.copy()
        z = xyz.copy()

        def numpy_assigment():
            z[:] = xyz

        print("normal method")
        print(timeit(normal, number=1000))
        print("use_memcpy")
        print(timeit(use_memcpy, number=1000))
        print("use_numpy")
        print(timeit(numpy_assigment, number=1000))

    @test_if_path_exists("./data/nogit/tip3p/")
    def test_2(self):
        print("tip3p data")
        print("Frame")
        mydir = "./data/nogit/tip3p/"
        traj = mdio.iterload(mydir + "/md.trj", mydir + "/tc5bwat.top")
        f0, f1 = traj[:2]
        xyz = f0.xyz + 1.0
        f0._fast_copy_from_xyz(xyz)
        aa_eq(f0.xyz, xyz)
        f1._fast_copy_from_frame(f0)
        aa_eq(f0.xyz, f1.xyz)

        frame = f0.copy()
        xyz = frame.xyz.copy()
        z = xyz.copy()

        def normal():
            frame.buffer2d[:] = xyz

        frame = f0.copy()
        xyz = frame.xyz.copy()
        z = xyz.copy()

        def use_memcpy():
            frame._fast_copy_from_xyz(xyz)

        frame = f0.copy()
        xyz = frame.xyz.copy()
        z = xyz.copy()

        def numpy_assigment():
            z[:] = xyz

        print("normal method")
        print(timeit(normal, number=1000))
        print("use_memcpy")
        print(timeit(use_memcpy, number=1000))
        print("use_numpy")
        print(timeit(numpy_assigment, number=1000))

        # RESULT: similiar copying times for all 3 methods
        # --> IO is bottle neck for all?


if __name__ == "__main__":
    unittest.main()
