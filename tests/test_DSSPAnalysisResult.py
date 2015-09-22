from __future__ import print_function
import unittest
import pytraj as pt
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        # just want to test if no error raised
        p = pt.load_pdb_rcsb("1l2y")
        d = pt.dssp(p, dtype='_dssp_class')
        d.average()
        d.to_ndarray()
        d.to_ndarray('string')
        d.to_ndarray('int')
        d.to_dict('int')
        d.to_dict('string')
        #print(d.residues)
        #print(np.vstack((d.residues, d.to_ndarray('string').T)))
        d.values_per_residue()
        d.values_per_frame()
        d.to_ndarray_per_frame()


if __name__ == "__main__":
    unittest.main()
