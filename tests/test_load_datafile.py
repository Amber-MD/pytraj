from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ds_saved = mdio.load_datafile("./data/dssp.Tc5b.dat")
        ds = traj.calc_dssp(dtype='dataset')
        print(sorted(ds.keys()))
        print(sorted(ds_saved.keys()))
        aa_eq(ds.filter("TYR:3").to_ndarray(),
              ds_saved.filter("TYR:3").to_ndarray())

        for i in range(2, 20):
            i_s = ":" + str(i)
            aa_eq(ds.filter(i_s).to_ndarray(),
                  ds_saved.filter(i_s).to_ndarray())

if __name__ == "__main__":
    unittest.main()
