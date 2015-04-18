from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir

test_dir = cpptraj_test_dir + "/Test_Comprehensive"

class Test(unittest.TestCase):
    @test_if_path_exists(cpptraj_test_dir)
    @test_if_having("pandas")
    def test_0(self):
        from pytraj import to_dataframe
        import os
        trajin_file = "./data/ptraj_comp.in"
        command = "cat %s" % trajin_file
        state = mdio.load_cpptraj_file(trajin_file)
        state.run()
        dslist = state.datasetlist
        print (dslist.get_legends())
        print (to_dataframe(dslist))


if __name__ == "__main__":
    unittest.main()
