from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = pt.search_hbonds(traj, )

        # +=
        d0 = dslist[0]
        arr0 = d0.copy().values
        d0.values += 1
        aa_eq(d0.data, arr0 + 1)

        # +
        arr0 = d0.copy().values
        d1 = d0 + 1.
        aa_eq(d1, arr0 + 1)

        # += for DatasetList
        arr0 = dslist[0].copy().values
        dslist[0].values += 2
        aa_eq(dslist[0].tolist(), arr0 + 2)

        # *=
        arr0 = dslist[0].copy().values
        dslist[0].values *= 2
        aa_eq(dslist[0].tolist(), arr0 * 2)

        # /=
        arr0 = dslist[0].to_ndarray().copy()
        # dslist[0] /= 2. # fail in PY3 because dslist[0] is DataSet_integer
        dslist[0] /= 2
        aa_eq(dslist[0].tolist(), arr0 / 2.)

        d = dslist[0] * 2 + 1

        d = dslist[0] * 2 + 1


if __name__ == "__main__":
    unittest.main()
