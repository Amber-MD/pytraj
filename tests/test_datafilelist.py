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
        from pytraj import DataFileList
        df = DataFileList()
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = traj.calc_dssp(dtype='dataset')

        for d in dslist:
            d.legend = d.legend.replace(":", "_")
            df.add_dataset("./output/testdf_" + d.legend + ".txt", d)
        df.write_all_datafiles()

        dslist.groupby("SER").write_all_datafiles()

if __name__ == "__main__":
    unittest.main()
