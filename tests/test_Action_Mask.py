import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        act = adict['mask']
        dslist = DataSetList()
        dflist = DataFileList()
        act('mask "(:5 <:3.0) & :WAT"', traj, dslist=dslist,
            dflist=dflist)
        print(dslist.size)
        # dflist.write_all_datafiles()

if __name__ == "__main__":
    unittest.main()
