from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import load_sample_data
from pytraj.utils import has_
from pytraj.testing import no_test


class Test(unittest.TestCase):
    @no_test
    def test_0(self):
        traj = load_sample_data()
        act = adict['matrix']
        dslist = DataSetList()

        # all atom distance matrix
        act("", traj, dslist=dslist)
        if has_("matplotlib"):
            from pytraj.plotting import plot_matrix
            outfit = plot_matrix(dslist[0])
            print(outfit)
            assert isinstance(outfit, tuple)
            for out in outfit:
                print(out)
        else:
            print("does not have matplotlib. skip test")


if __name__ == "__main__":
    unittest.main()
