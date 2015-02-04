import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.cast_dataset import cast_dataset

farray = TrajReadOnly(top=Topology("./data/Tc5b.top"), 
                    filename='data/md1_prod.Tc5b.x', 
                    )

class TestDistance(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        dslist2 = DataSetList()
        act = allactions.Action_Distance()
        act_radgyr = allactions.Action_Radgyr()

        act.master(command="distance :1@CA :2@CA", 
                   traj=farray,
                   current_top=farray.top, dslist=dslist)

        act_radgyr.master(command="radgyr @CA",
                   traj=farray,
                   current_top=farray.top, dslist=dslist2)
       
        d1 = cast_dataset(dslist[0], dtype="general")
        d2 = cast_dataset(dslist2[0], dtype="general")
        print(d1.data[:10])
        print(d2.data[:10])

if __name__ == "__main__":
    unittest.main()
