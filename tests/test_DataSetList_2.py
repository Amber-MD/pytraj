import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets import cast_dataset

class Test(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ds = dslist.add_set_aspect('INTEGER', 'frame_idx', 'frame_idx')
        ds0 = cast_dataset(ds, 'INTEGER')
        print (dir(ds0))

        for i in range(100):
            ds0.add(i, i)

        print (dslist[0].size)
        print (dslist[0][:])
        print (dslist[0].aspect)
        dslist2 = DataSetList()
        act = adict['distance']
        act(":2@CA :10@CA", traj, dslist=dslist2)
        print (dslist2[0].aspect)

        print ("get_xyz methods")
        print (dslist2.get_legends())
        print (dslist2.get_aspects())
        print (dslist2.get_scalar_modes())
        print (dslist2.get_scalar_types())
        print (dslist2.get_dtypes())

        for d0 in dslist2:
            print (d0.legend)
            print (d0.aspect)
            print (d0.scalar_mode)
            print (d0.scalar_type)


if __name__ == "__main__":
    unittest.main()
