import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.cast_dataset import cast_dataset

class Test(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        ds = dslist.add_set_aspect('INTEGER', 'frame_idx', 'frame_idx')
        ds0 = cast_dataset(ds, 'INTEGER')
        print (dir(ds0))

        for i in range(100):
            ds0.add(i, i)

        print (dslist[0].size)
        print (dslist[0][:])

if __name__ == "__main__":
    unittest.main()
