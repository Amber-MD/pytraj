import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.cast_dataset import cast_dataset

farray = TrajReadOnly(top=Topology("./data/Tc5b.top"), 
                    filename='data/md1_prod.Tc5b.x', 
                    )

class TestRadgyr(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        act = allactions.Action_Radgyr()
        act.read_input("radgyr @CA", farray.top, dslist=dslist)
        act.process(farray.top)
       
        for i, frame in enumerate(farray):
            act.do_action(i, frame)

        d1 = cast_dataset(dslist[0], dtype="general")
        print(d1.data[:10])
        print(dir(d1))
        #print d1.size
        #print cast_dataset.__doc__

if __name__ == "__main__":
    unittest.main()
