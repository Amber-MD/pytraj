from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.misc import to_amber_mask
        assert to_amber_mask('ASP_16@OD1-ARG_18@N-H') == ":16@OD1 :18@N"
        saved_mask_list = [":16@OD1 :18@N", ":16@OD1 :18@N"]
        assert to_amber_mask(
            ['ASP_16@OD1-ARG_18@N-H', 'ASP_16@OD1-ARG_18@N-H'
             ]) == saved_mask_list
        print(
            to_amber_mask(['ASP_16@OD1-ARG_18@N-H', 'ASP_16@OD1-ARG_18@N-H']))


if __name__ == "__main__":
    unittest.main()
