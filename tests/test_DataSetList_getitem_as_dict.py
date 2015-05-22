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

        dslist2 = DataSetList()
        act = adict['distance']
        act(":2@CA :10@CA", traj, dslist=dslist2)
        print (dslist2.get_legends())
        assert (dslist2['Dis_00000'].size == 10)

        # getting a list of legends
        dslist3 = traj.search_hbonds()

        legends = [dslist3[5].legend, dslist3[4].legend]
        print (dslist3[legends])

if __name__ == "__main__":
    unittest.main()
