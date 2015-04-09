"""This file is dedicated for repeating the same command, actions ...
to make sure there is no segmentation fault, no memory error"""
from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having

class Test(unittest.TestCase):

    #@no_test
    def test_0(self):
        print ("repeat `calculate`")
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        from pytraj import calculate

        d0 = calculate("distance", ":2@CA :10@CB", traj)
        print (d0.size)
        d1 = calculate("distance", ":3@CA :10@CB", traj)
        print (d1.size)
        #assert d0.size == d1.size == traj.n_frames

    #@no_test
    def test_1(self):
        print ("repeat `calc_`")
        str_traj_top = ("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj = mdio.load(*str_traj_top)
        from pytraj.common_actions import calc_distance

        d0 = calc_distance(":2@CA :10@CA", traj)
        d1 = calc_distance(":3@CA :10@CB", traj)
        d2 = calc_distance(":3@CA :10@CB", traj)
        d3 = calc_distance(":3@CA :10@CB", traj)
        d4 = calc_distance(":3@CA :10@CB", *str_traj_top)

        assert d0.size == d1.size == traj.n_frames 
        assert d0.size == d2.size == d3.size
        assert d0.size == d4.size

if __name__ == "__main__":
    unittest.main()
