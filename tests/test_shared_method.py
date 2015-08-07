from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj._shared_methods import _frame_iter


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.trajs.Trajin_Single import Trajin_Single
        traj0 = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj = Trajin_Single(traj0.filename, traj0.top)

        print("TrajectoryIterator")
        print(traj)
        for frame in _frame_iter(traj, 1, 8, 2, '@CA'):
            print(frame)

        print("Trajectory")
        farray = traj[:]
        print(farray)
        for frame in _frame_iter(farray, 1, 8, 2, '@CA'):
            print(frame)

        print("CRD dataset")
        from pytraj.datasets import DataSet_Coords_CRD
        crd = DataSet_Coords_CRD()
        crd.load(farray, farray.top)
        print(crd)
        for frame in _frame_iter(crd, 1, 8, 2, '@CA'):
            print(frame)

        print("TRAJ dataset")
        from pytraj.datasets import DataSet_Coords_TRJ
        coords_traj = DataSet_Coords_TRJ()
        coords_traj.add_trajin(traj)
        print(coords_traj)
        for frame in _frame_iter(coords_traj, 1, 8, 2, '@CA'):
            print(frame)


if __name__ == "__main__":
    unittest.main()
