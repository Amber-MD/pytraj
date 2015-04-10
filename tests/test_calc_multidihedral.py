from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having

class Test(unittest.TestCase):
    @test_if_having("pandas")
    def test_0(self):
        from pytraj.dataframe import to_dataframe
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        from pytraj.common_actions import calc_multidihedral
        command =  "resrange 2-19 phi psi"
        d0 = calc_multidihedral(command, traj)
        d1 = calc_multidihedral(command, traj)
        assert isinstance(d0, dict) == True
        assert (len(d0.keys()) == len(d1.keys()))
        print (to_dataframe(d1))
        print (to_dataframe(d0))

        d3 = calc_multidihedral("", traj)
        print (to_dataframe(d3))

if __name__ == "__main__":
    unittest.main()
