from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        mask = pt.utils.convert.array2d_to_cpptraj_maskgroup([[0, 3], [5, 6, 8]
                                                              ])
        assert mask == '@1,4 @6,7,9'


if __name__ == "__main__":
    unittest.main()
