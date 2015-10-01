from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca


class TestGrid(unittest.TestCase):
    def test_0(self):
        from numpy.testing import assert_almost_equal as a_equal
        from pytraj.math import Grid
        import numpy as np
        nx = ny = nz = 3
        g = Grid(nx, ny, nz)
        assert g.size == nx ** 3
        assert g.nx == g.ny == g.nz == nx

        value = 1000.
        g[0, 0, 0] = value
        assert g[0, 0, 0] == value
        assert g._element(0, 0, 0) == value

        np_arr = g.to_ndarray()
        a_list = g.tolist()
        a_equal(np_arr, a_list)


class TestGridAction(unittest.TestCase):
    def test_0(self):
        traj = pt.load_sample_data("tz2")[:]
        traj.autoimage()
        traj.rmsfit(mask=':1-13')
        d = pyca.calc_grid(traj, " 20 0.5 20 0.5 20 0.5 :WAT@O")

        d = pyca.calc_grid(
            traj(), " 20 0.5 20 0.5 20 0.5 :WAT@O",
            top=traj.top)


if __name__ == "__main__":
    unittest.main()
