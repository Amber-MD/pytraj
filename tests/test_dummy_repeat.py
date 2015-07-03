"""This file is dedicated for repeating the same command, actions ...
to make sure there is no segmentation fault, no memory error"""
from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having
import pytraj.common_actions as pyca
from pytraj.compat import range


class Test(unittest.TestCase):

    def test_0(self):
        print("repeat `calculate`")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        from pytraj import calculate

        d0 = calculate("distance", traj, ":2@CA :10@CB")
        print(d0.size)
        d1 = calculate("distance", traj, ":3@CA :10@CB")
        print(d1.size)
        assert d0[0].size == d1[0].size == traj.n_frames

    def test_1(self):
        print("repeat `calc_`")
        str_traj_top = ("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj = mdio.iterload(*str_traj_top)
        from pytraj.common_actions import calc_distance

        d0 = calc_distance(traj, ":2@CA :10@CA")
        print(d0.size, d0)
        d1 = calc_distance(traj, ":3@CA :10@CB")
        print(d1.size)
        d2 = calc_distance(traj, ":3@CA :10@CB")
        d3 = calc_distance(traj, ":3@CA :10@CB")

        assert d0.size == d1.size == traj.n_frames
        assert d0.size == d2.size == d3.size

    def test_2(self):
        # calc_dssp
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # repeat 20 times
        for i in range(20):
            d0 = pyca.calc_dssp(traj[:], dtype='dataset')
            assert d0.size == 28
            # WARNING: don't assign `d0 = d0.filter("your_mask")`
            d1 = d0.filter("TRP:6")

        print(d1.keys())

    def test_3(self):
        print("repeat to_dict")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d0 = traj.search_hbonds(dtype='dataset')
        ddcit = d0.to_dict()

        for i in range(100):
            d = pyca.search_hbonds(traj, dtype='dataset').to_dict()
            assert sorted(ddcit) == sorted(d)

        assert (ddcit.keys().__len__() == 20)

if __name__ == "__main__":
    unittest.main()
