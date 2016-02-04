from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        from itertools import product
        import numpy as np
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        d0 = pt.calc_mindist(traj, "@CA @CB")
        i0 = top("@CA").indices
        i1 = top("@CB").indices
        combinations = np.array(list(product(i0, i1)))
        d1 = pt.distance(traj, combinations).T

        min_list = []
        for arr0 in d1:
            min_list.append(np.min(arr0))
        aa_eq(d0, min_list)

        d0 = pt.mindist(traj, [[0, 1], [2, 3]])
        d1 = pt.mindist(traj, '@1,2 @3,4')
        aa_eq(d0, d1)


if __name__ == "__main__":
    unittest.main()
