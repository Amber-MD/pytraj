from __future__ import print_function
import pytraj as pt
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir
from pytraj.datasetlist import stack

traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")


class Test(unittest.TestCase):

    def test_0(self):
        import numpy as np
        from pytraj.datasetlist import stack as stack
        _ds1 = pt.calc_dssp(traj[:5], dtype='dataset')
        ds1 = _ds1.grep('LYS')
        _ds2 = pt.calc_dssp(traj[5:], dtype='dataset')
        ds2 = _ds2.grep('LYS')

        dstack = stack((ds1, ds2))
        _d12 = pt.calc_dssp(traj, dtype='dataset')
        d12 = _d12.grep("LYS")

        dstack_dict = dstack.to_dict()
        d12_dict = d12.to_dict()
        assert sorted(dstack_dict.keys()) == sorted(d12_dict)

        for key in dstack_dict.keys():
            arr0 = dstack_dict[key]
            arr1 = d12_dict[key]
            if np.any(arr0 == arr1) == False:
                pass

        arr1 = ds1.to_ndarray()
        arr2 = ds2.to_ndarray()
        arrstack = dstack.to_ndarray()
        arr12 = d12.to_ndarray()

        aa_eq(arrstack.flatten(), arr12.flatten())

    def test_1(self):
        dslist0 = pt.calc_phi(traj)
        dslist1 = pt.calc_psi(traj)
        dslist2 = pt.search_hbonds(traj)
        self.assertRaises(KeyError, lambda: stack((dslist0, dslist1)))
        self.assertRaises(TypeError, lambda: stack((dslist0, dslist2)))

        stack((dslist0 for _ in range(3)))


if __name__ == "__main__":
    unittest.main()
