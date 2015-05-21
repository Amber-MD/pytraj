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
        dslist = traj.search_hbonds()

        # +=
        d0 = dslist[0]
        arr0 = d0.to_ndarray().copy()
        d0 += 1
        aa_eq(d0.data, arr0 + 1)

        # +
        arr0 = d0.to_ndarray().copy()
        d1 = d0 + 1.
        print (d1.tolist())
        aa_eq(d1.data, arr0 + 1)

        # += for DataSetList
        arr0 = dslist[0].to_ndarray().copy()
        dslist[0] += 2
        aa_eq(dslist[0].tolist(), arr0 + 2)

        # *=
        arr0 = dslist[0].to_ndarray().copy()
        dslist[0] *= 2
        aa_eq(dslist[0].tolist(), arr0 * 2)

        # /=
        arr0 = dslist[0].to_ndarray().copy()
        dslist[0] /= 2.
        print (dslist[0].tolist())
        print (arr0/2)
        print (arr0/2.)
        aa_eq(dslist[0].tolist(), arr0 / 2.)

        d = dslist[0] * 2 + 1
        print (d.tolist())

        d = dslist[0] * 2 + 1
        print (d.tolist())

if __name__ == "__main__":
    unittest.main()
