from __future__ import print_function
import unittest


class Test(unittest.TestCase):
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
        #print(g._element(0, 0, 0))

        #print(g._element(0, 0, 5))
        np_arr = g.to_ndarray()
        a_list = g.tolist()
        a_equal(np_arr, a_list)
        #print(np.asarray(g.data))


if __name__ == "__main__":
    unittest.main()
