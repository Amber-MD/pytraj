from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        # just want to test if no error raised
        p = pt.load_pdb_rcsb("1l2y")
        d = p.calc_dssp(dtype='_dssp_class')
        d.average()
        d.to_ndarray()
        d.to_ndarray('string')
        d.to_ndarray('int')
        d.to_dict('int')
        d.to_dict('string')


if __name__ == "__main__":
    unittest.main()
