import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.cpptraj_dict import DataTypeDict as ddict
from pytraj.cpptraj_dict import scalarDict, scalarModeDict

class Test(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        dflist = DataFileList()
        dist = dslist.add_set('double', "myname", "dis_")
        print (scalarModeDict.keys())
        dist.set_scalar('m_distance')

if __name__ == "__main__":
    unittest.main()
