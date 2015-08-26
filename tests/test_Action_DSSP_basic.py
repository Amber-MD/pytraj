import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.datasets import cast_dataset
from pytraj import adict
from pytraj.dssp_analysis import to_string_ss
from pytraj.core.DataFileList import DataFileList

traj = TrajectoryIterator(top=Topology("./data/DPDP.parm7"),
                          filename='./data/DPDP.nc', )


class TestRadgyr(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        act = adict['dssp']
        dflist = DataFileList()
        act(":10-22 out ./output/_test_dssp_DPDP.dat", traj,
            dslist=dslist,
            dflist=dflist)
        #print(dslist.size)

        d0 = dslist[-1]
        arr0 = dslist.get_dataset(dtype="integer")
        #print(arr0)

        from pytraj.utils import _import
        has_np, np = _import('numpy')

        if has_np:
            arr0_np = np.asarray(arr0)
            #print(arr0_np)

        for d0 in dslist:
            #print(d0.name)

        #print(dslist.get_legends())
        #print(dslist['LYS:17'])
        #print(dslist['LYS:17'][:])
        #print(dslist['DSSP_00000[Anti]'])
        #print(dslist['DSSP_00000[Anti]'][:])


if __name__ == "__main__":
    unittest.main()
