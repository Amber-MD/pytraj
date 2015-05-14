from __future__ import print_function
import unittest
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = traj.search_hbonds()
        mystr = dslist[0].__str__()

        expected_str0 = """<pytraj.datasets.DataSet_integer: size=10, name=HB_00000,"""
        expected_str1 = """ legend=HB_00000[UU], aspect=UU, dtype=integer, data_format= %12i>"""

        mystr = mystr.replace("\n", "")
        expected_str = expected_str0 + expected_str1

        assert mystr == expected_str

if __name__ == "__main__":
    unittest.main()
