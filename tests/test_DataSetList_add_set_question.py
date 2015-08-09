from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets.DataSet_Coords_TRJ import DataSet_Coords_TRJ
from pytraj.datasets.DatasetMatrixDouble import DatasetMatrixDouble


class Test(unittest.TestCase):
    # FIXME: need to add dataset's name to add to DSL

    def test_0(self):
        dset_traj = DataSet_Coords_TRJ()
        dslist = DataSetList()

        # wrapper of "AddSet(DataSet*)
        dslist.add_existing_set(dset_traj)
        print(dslist.size)  # = 0, but I expected "=1"

    def test_1(self):
        dset_traj = DatasetMatrixDouble()
        dslist = DataSetList()

        # wrapper of "AddSet(DataSet*)
        dslist.add_existing_set(dset_traj)
        print(dslist.size)  # = 0, but I expected "=1"

    def test_2(self):
        dslist = DataSetList()

        # wrapper of "AddSet(DataType, name, default_name)"
        dslist.add_set("coords", "name", "funny_name")
        print(dslist[0])
        print(dslist.size)  # = 1 (is what I expected)


if __name__ == "__main__":
    unittest.main()
