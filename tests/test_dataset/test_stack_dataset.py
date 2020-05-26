from __future__ import print_function
import numpy as np
import pytraj as pt
from utils import fn
import unittest

from pytraj import io as mdio
from pytraj.testing import aa_eq
from pytraj.datasets.datasetlist import stack
from pytraj.datasets.datasetlist import stack as stack
import pytest

traj = mdio.iterload(fn('Tc5b.x'), fn('Tc5b.top'))


class Test(unittest.TestCase):
    def test_0(self):
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
        ds2.to_ndarray()
        arrstack = dstack.to_ndarray()
        arr12 = d12.to_ndarray()

        aa_eq(arrstack.flatten(), arr12.flatten())

    def test_1(self):
        dslist0 = pt.calc_phi(traj)
        dslist1 = pt.calc_psi(traj)
        dslist2 = pt.search_hbonds(traj)
        with pytest.raises(KeyError):
            stack((dslist0, dslist1))
        with pytest.raises(TypeError):
            stack((dslist0, dslist2))

        stack((dslist0 for _ in range(3)))


if __name__ == "__main__":
    unittest.main()
