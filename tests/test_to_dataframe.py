from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.dataframe import to_dataframe, has_pandas
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['multidihedral']
        dslist = DataSetList()
        act("phi", traj, dslist=dslist)
        #print (dslist)
        
        if has_pandas:
            print ("has_pandas")
            dframe = to_dataframe(dslist)
            print (dframe)
        else:
            print ("does not have pandas installed")
            print ("skip")

if __name__ == "__main__":
    unittest.main()
