from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
import numpy as np
from numpy.testing import assert_almost_equal as aa_eq

"""enhance DataSetList by using numpy"""

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ds = traj.search_hbonds()
        arr0 = ds.to_ndarray()
        aa_eq(arr0.mean(1), ds.mean())
        aa_eq(np.median(arr0, 1), ds.median())
        aa_eq(arr0.sum(1), ds.sum())
        aa_eq(arr0.min(1), ds.min())
        aa_eq(arr0.max(1), ds.max())
        aa_eq(arr0.std(1), ds.std())
        ds.apply(lambda x: 2*x)
        newarr = ds.to_ndarray()
        aa_eq(2*arr0, newarr)

if __name__ == "__main__":
    unittest.main()
