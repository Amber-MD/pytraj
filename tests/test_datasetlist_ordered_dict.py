from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = traj.calc_multidihedral()
        assert np.all(
            pt.tools.dict_to_ndarray(dslist.to_dict(ordered_dict=True)) ==
            dslist.to_ndarray())

        self.assertRaises(NotImplementedError,
                          lambda: pt.tools.dict_to_ndarray(dslist.to_dict()))


if __name__ == "__main__":
    unittest.main()
