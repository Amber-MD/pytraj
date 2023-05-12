from __future__ import print_function
import unittest
import pytraj as pt
import numpy as np
from utils import fn


class Test(unittest.TestCase):
    def test_0(self):
        mask = pt.utils.convert.array2d_to_cpptraj_maskgroup(np.array([[0, 3],
                                                              [5, 6, 8]], dtype=object))
        assert mask == '@1,4 @6,7,9'


if __name__ == "__main__":
    unittest.main()
