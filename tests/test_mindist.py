from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        from itertools import product
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        d0 = pyca.calc_mindist(traj, "@CA @CB")
        i0 = top("@CA").indices
        i1 = top("@CB").indices
        combinations = np.array(list(product(i0, i1)))
        d1 = traj.calc_distance(combinations)
        print (d1.shape)

        min_list = []
        for arr0 in d1:
            min_list.append(np.min(arr0))
        print (min_list)
        print (d0.tolist())
        aa_eq(d0.tolist(), min_list)

if __name__ == "__main__":
    unittest.main()
