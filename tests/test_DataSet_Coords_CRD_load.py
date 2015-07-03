from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets.DataSet_Coords_CRD import DataSet_Coords_CRD


class Test(unittest.TestCase):

    def test_0(self):
        print("load TrajectoryIterator")
        coords = DataSet_Coords_CRD()
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        coords.top = traj.top
        coords.load(traj)
        assert coords.size == traj.size

    def test_1(self):
        print("load Trajectory")
        coords = DataSet_Coords_CRD()
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        coords.top = traj.top
        coords.load(farray)
        assert coords.size == traj.size

    def test_2(self):
        print("load string (filenames)")
        coords = DataSet_Coords_CRD()
        coords.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        assert coords.size == traj.size

    def test_3(self):
        print("load frame_iter")
        coords = DataSet_Coords_CRD()
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        coords.load(traj(2, 8, 2), traj.top)
        assert coords.size == 3

if __name__ == "__main__":
    unittest.main()
