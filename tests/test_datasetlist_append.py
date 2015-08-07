from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        # test copy
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = traj.calc_multidihedral("psi", dtype='dataset')
        self.assertRaises(KeyError, lambda: dslist.append(dslist[0]))

        d0 = dslist[0].copy()
        d0.key = 'new_key0'
        dslist.append(d0, copy=True)
        assert d0 not in dslist

        d1 = dslist[0].copy()
        d1.key = 'new_key1'
        dslist.append(d1, copy=False)
        assert d1 in dslist

        # always make a copy with constructor
        dslist2 = dslist.__class__(dslist)

        for d0 in dslist:
            assert d0 not in dslist2

        # test copy=False for append
        dslist3 = dslist.__class__()
        for d0 in dslist:
            dslist3.append(d0, copy=False)
            assert d0 in dslist3

        # test grep copy=False
        dslist4 = dslist.grep("psi", copy=False)
        for d0 in dslist4:
            assert d0 in dslist

        # test grep copy=True
        dslist5 = dslist.grep("psi", copy=True)
        for d0 in dslist5:
            assert d0 not in dslist


if __name__ == "__main__":
    unittest.main()
