import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.cast_dataset import cast_dataset
from pytraj import adict 
from pytraj.DataFileList import DataFileList

farray = TrajReadOnly(top=Topology("./data/DPDP.parm7"), 
                    filename='./data/DPDP.nc', 
                    )

class TestRadgyr(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        act = adict['dssp']
        dflist = DataFileList()
        act.read_input(":2-15 out ./output/_test_dssp_DPDP.dat", farray.top, dslist=dslist, dflist=dflist)
        act.process(farray.top)
       
        for i, frame in enumerate(farray):
            act.do_action(i, frame)

        print (dslist.size)
        dflist.write_all_datafiles()

        print (dslist.size)
        for dset in dslist:
            print (dset.dtype)
            print (dset[:])

        # TODO : intepret the output (not understand what they mean)

if __name__ == "__main__":
    unittest.main()
