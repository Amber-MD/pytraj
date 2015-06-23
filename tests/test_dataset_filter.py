from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d0 = traj.calc_multidihedral(dtype='dataset')[0]
        import numpy as np
        arr0 = np.array([ 154.85298804,  130.30039598,  153.107096  ])
        print (d0)
        aa_eq(d0.filter(lambda x : 130. < x < 155.), arr0)


if __name__ == "__main__":
    unittest.main()
