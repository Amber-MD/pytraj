import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        # TODO : add assert
        # TODO : correct casting data type
        # make dictionary? (sound good)
        traj = mdio.load("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        act = adict['hbond']
        dslist = DataSetList()
        act(":1-13 solventacceptor :WAT@O solventdonor :WAT series", 
            traj, dslist=dslist)
        act.print_output()
        print ('dslist size = ', dslist.size)

        for d0 in dslist:
            if d0.dtype == 'integer':
                print (d0[:])
            print (d0.name, d0)

        for i in range(d0.size):
            print (i)
        act.help()

        d3 = dslist[3]
        print (d3)
        print ([d3[i] for i in range(d3.size)])
        print (dslist.get_dtypes())
        print (dslist.get_aspects())
        print (dslist.get_scalar_modes())
        print (dslist.get_scalar_types())
        print (dslist.get_legends())

if __name__ == "__main__":
    unittest.main()
