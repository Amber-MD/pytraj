from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.datasets.DataSetList import DataSetList as CpptrajDSL
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = CpptrajDSL()
        pt.adict['multidihedral']("", traj, top=traj.top, dslist=dslist)
        da = pt.array.DataArray(dslist[0])
        print(da)

        # test copy
        assert da.copy().values is not da.values
        assert da.shallow_copy().values is da.values

if __name__ == "__main__":
    unittest.main()
