import unittest
import numpy as np
import pytraj as pt
from pytraj import TrajectoryIterator
from pytraj.base import DataSetList, DataFileList
from numpy.testing import assert_allclose
from pytraj import io as mdio
from pytraj import allactions
from pytraj.datasets import cast_dataset
from pytraj import adict
from pytraj.core.DataFileList import DataFileList
from pytraj.common_actions import calc_dssp
from pytraj.testing import aa_eq

farray = TrajectoryIterator(top=pt.load_topology("./data/DPDP.parm7"),
                            filename='./data/DPDP.nc', )


class TestRadgyr(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        act = adict['dssp']
        dflist = DataFileList()
        act.read_input(":10-22 out ./output/_test_dssp_DPDP.dat", farray.top,
                       dslist=dslist,
                       dflist=dflist)
        act.process(farray.top)

        for i, frame in enumerate(farray):
            act.do_action(frame)

        dflist.write_all_datafiles()

        arr1 = dslist.get_dataset(dtype='float')
        arr0 = dslist.get_dataset(dtype='integer')

    def test_1(self):
        dslist = DataSetList()
        dflist = DataFileList()
        adict['dssp'](":10-22 out ./output/_test_dssp_DPDP.dat",
                      current_frame=farray,
                      top=farray.top,
                      dslist=dslist,
                      dflist=dflist)
        # Secondary structure for each residue in mask for 100 frames

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


if __name__ == "__main__":
    unittest.main()
