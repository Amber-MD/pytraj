import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        # TODO : add assert
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['hbond']
        dslist = DataSetList()
        act("series", traj, dslist=dslist)
        print ('dslist size = ', dslist.size)
        for d0 in dslist:
            print (d0[:])
            print (d0.name)

        act.help()

if __name__ == "__main__":
    unittest.main()
