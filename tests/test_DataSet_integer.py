import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.datasets.DataSet_integer import DataSet_integer
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ds = DataSet_integer()
        print (ds.dtype)

        for i in range(100):
            ds.add_element(i)

        assert ds.size == 100
        print (ds[:])
        print (dir(ds))

if __name__ == "__main__":
    unittest.main()
