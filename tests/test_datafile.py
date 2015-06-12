from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        # basic (long)
        from pytraj.core.DataFileList import DataFileList
        from pytraj import ArgList
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ds = pyca.calc_molsurf(traj)
        print (ds.legend)

        df = DataFileList()
        d0 = df.add_datafile("output/test_datafile.txt", ArgList())
        d0.add_dataset(ds)
        d0.write_data()

        # from DataSet
        ds.write_to_cpptraj_format("./output/test_datafile_cpptraj.txt")

if __name__ == "__main__":
    unittest.main()
