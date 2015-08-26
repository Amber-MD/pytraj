from __future__ import print_function
import unittest; import pytraj as pt
from pytraj import io as mdio

expected_str_hbonds = """<pytraj.datasets.DataSet_integer: size=10, name=HB_00000, legend=HB_00000[UU], aspect=UU, dtype=integer, format= %12i>
values: [0 3 3 2 2 3 4 4 2 2]"""


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = pt.search_hbonds(traj,)
        mystr = dslist[0].__str__()
        #print(mystr)


if __name__ == "__main__":
    unittest.main()
