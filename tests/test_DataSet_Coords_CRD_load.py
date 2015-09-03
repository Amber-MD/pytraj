from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets.DataSet_Coords_CRD import DataSet_Coords_CRD


class Test(unittest.TestCase):
    def test_0(self):
        coords = DataSet_Coords_CRD()
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        coords.top = traj.top
        coords.load(traj)
        assert coords.n_frames == traj.n_frames

    def test_1(self):
        coords = DataSet_Coords_CRD()
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        coords.top = traj.top
        coords.load(farray)
        assert coords.n_frames == traj.n_frames

    def test_2(self):
        coords = DataSet_Coords_CRD()
        coords.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        assert coords.n_frames == traj.n_frames

    def test_3(self):
        coords = DataSet_Coords_CRD()
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        coords.load(traj(2, 8, 2), traj.top)
        assert coords.n_frames == 3


if __name__ == "__main__":
    unittest.main()
