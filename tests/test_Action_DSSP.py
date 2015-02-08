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
        act.read_input(":10-22 out ./output/_test_dssp_DPDP.dat", farray.top, dslist=dslist, dflist=dflist)
        act.process(farray.top)
       
        for i, frame in enumerate(farray):
            act.do_action(i, frame)

        print (dslist.size)
        dflist.write_all_datafiles()

        print (dslist.size)
        for dset in dslist:
            print (dset.dtype)

        # TODO : intepret the output (not understand what they mean)
        arr1 = dslist.get_dataset(dtype='float')
        print (arr1.shape)

        arr0 = dslist.get_dataset(dtype='integer')
        print (arr0.shape)

        print (arr0[0].__len__())
        print (arr0[0])

if __name__ == "__main__":
    unittest.main()
