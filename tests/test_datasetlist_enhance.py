from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
import numpy as np
from numpy.testing import assert_almost_equal as aa_eq
from pytraj.datasets.DataSetList import DataSetList

"""enhance DataSetList by using numpy"""


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ds = traj.search_hbonds()
        arr0 = ds.to_ndarray()
        aa_eq(arr0.mean(1), ds.mean()[1])
        aa_eq(arr0.sum(1), ds.sum()[1])
        aa_eq(arr0.min(1), ds.min()[1])
        aa_eq(arr0.max(1), ds.max()[1])
        aa_eq(arr0.std(1), ds.std()[1])

        aa_eq(arr0.mean(1),ds.mean()[1])
        aa_eq(arr0.sum(1), ds.sum()[1])
        aa_eq(arr0.min(1), ds.min()[1])
        aa_eq(arr0.max(1), ds.max()[1])
        aa_eq(arr0.std(1), ds.std()[1])

        ds.apply(lambda x: 2 * x)
        newarr = ds.to_ndarray()
        aa_eq(2 * arr0, newarr)

    def test_1(self):
        print("test __array__")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ds = traj.search_hbonds()
        ds0 = ds[0]
        print(ds0)
        assert np.mean(ds0.values) == ds0.avg()
        assert np.sum(ds0) == np.sum(ds0.data)

        print("test split")
        print(ds0)
        print(ds0.split(3))

if __name__ == "__main__":
    unittest.main()
