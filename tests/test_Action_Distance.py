import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.datasets import cast_dataset

farray = TrajectoryIterator(top=Topology("./data/Tc5b.top"),
                            filename='data/md1_prod.Tc5b.x',
                            )


class TestDistance(unittest.TestCase):

    def test_0(self):
        dslist = DataSetList()
        dslist2 = DataSetList()
        act = allactions.Action_Distance()
        act_radgyr = allactions.Action_Radgyr()

        act(command="distance :1@CA :2@CA",
            current_frame=farray,
            top=farray.top, dslist=dslist)

        act_radgyr(command="radgyr @CA",
                   current_frame=farray,
                   top=farray.top, dslist=dslist2)

        d1 = cast_dataset(dslist[0], dtype="general")
        d2 = cast_dataset(dslist2[0], dtype="general")
        print(d1.data[:10])
        print(d2.data[:10])

if __name__ == "__main__":
    unittest.main()
