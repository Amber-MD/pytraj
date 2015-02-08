import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.cast_dataset import cast_dataset
from pytraj import adict 
from pytraj.DataFileList import DataFileList

farray = TrajReadOnly(top=Topology("./data/Tc5b.top"), 
                    filename='data/md1_prod.Tc5b.x', 
                    )

class TestRadgyr(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        act = adict['dssp']
        dflist = DataFileList()
        act.read_input("out test_dssp.out", farray.top, dslist=dslist, dflist=dflist)
        act.process(farray.top)
       
        for i, frame in enumerate(farray):
            act.do_action(i, frame)

        print (dslist.size)
        dflist.write_all_datafiles()
        print (dslist[0])

if __name__ == "__main__":
    unittest.main()
