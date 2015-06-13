import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.datasets.DatasetInteger import DatasetInteger
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ds = DatasetInteger()
        print (ds.dtype)

        for i in range(100):
            ds.append(i)

        assert ds.size == 100
        print (ds[:])
        print (dir(ds))
        ds2 = DatasetInteger()

        ds2.append(ds) # copy from ds
        aa_eq(ds2.data, ds.data)

if __name__ == "__main__":
    unittest.main()
