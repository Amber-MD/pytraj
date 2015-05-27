import unittest
import numpy as np
from numpy.testing import assert_allclose
from pytraj import io as mdio
from pytraj.base import *
from pytraj import allactions
from pytraj.datasets import cast_dataset
from pytraj import adict 
from pytraj.DataFileList import DataFileList
from pytraj.common_actions import calc_dssp
from pytraj.testing import aa_eq

farray = TrajectoryIterator(top=Topology("./data/DPDP.parm7"), 
                    filename='./data/DPDP.nc', 
                    )

class TestRadgyr(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        act = adict['dssp']
        dflist = DataFileList()
        act.read_input(":10-22 out ./output/_test_dssp_DPDP.dat", 
                      farray.top, dslist=dslist, dflist=dflist)
        act.process(farray.top)
       
        for i, frame in enumerate(farray):
            act.do_action(frame)

        print (dslist.size)
        dflist.write_all_datafiles()

        print (dslist.size)
        for dset in dslist:
            print (dset.dtype)

        arr1 = dslist.get_dataset(dtype='float')
        arr0 = dslist.get_dataset(dtype='integer')
        print (arr0[0].__len__())
        print (arr0[0])

    def test_1(self):
        dslist = DataSetList()
        dflist = DataFileList()
        adict['dssp'](":10-22 out ./output/_test_dssp_DPDP.dat", 
            current_frame=farray, top=farray.top, 
            dslist=dslist, dflist=dflist)
        print (dslist.size)
        arr0 = dslist.get_dataset(dtype="integer")

        # Secondary structure for each residue in mask for 100 frames

    def test_2(self):
        def calc_dssp(command="", traj=None):
            dslist = DataSetList()
            adict['dssp'](command, 
                          current_frame=traj, top=traj.top, 
                          dslist=dslist)
            return dslist.get_dataset(dtype="integer")
        arr0 = calc_dssp(":10-22", farray[:2])

    def test_3(self):
        arr0 = calc_dssp(farray[:3], ":10-22", dtype='int')
        arr1 = calc_dssp(farray[:3], ":10-22", dtype='str')
        print (arr0)
        print (arr1)

    def test_4(self):
        # add assert 
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        py_d = calc_dssp(traj, "*", dtype='dataset')
        # load cpptraj output
        cpp_d = mdio.load_datafile("./data/dssp.Tc5b.dat")
        for key in cpp_d.keys():
            aa_eq(py_d[key].to_ndarray(), cpp_d[key].to_ndarray())

    def test_5(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = calc_dssp(traj, "*", dtype='dataset')
        print (dslist)
        print (dslist.get_legends())

if __name__ == "__main__":
    unittest.main()
