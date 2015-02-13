import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from glob import glob
from pytraj import adict

class Test(unittest.TestCase):
    def test_0(self):
        datadir = "./data/Test_RemdTraj/"
        dfiles = sorted(glob(datadir + "rem.nc.*"))
        print (dfiles)
        top = Topology(datadir + "ala2.99sb.mbondi2.parm7")

        trajlist = [mdio.load(traj, top) for traj in dfiles]
        trajlist2 = [mdio.load(traj, top)(1, 8, 2) for traj in dfiles]
        trajlist3 = [mdio.load(traj, top)(start=5) for traj in dfiles]

        print (trajlist)
        print (trajlist[0])
        dslist = DataSetList()
        adict['radgyr']("@CA", trajlist, top, dslist=dslist)
        print (dslist.size)
        print (dslist[0].size)

        d0 = adict['radgyr']("@C*", trajlist2, top, quick_get=True)
        print (d0.size)

        d0 = adict['radgyr']("@C*", trajlist3, top, quick_get=True)
        print (d0.size)

        top2 = top.strip_atoms("!@CA", copy=True)
        d0 = adict['radgyr']("@C*", [traj['@CA :frame'] for traj in trajlist], 
                             top2, quick_get=True)
        print (d0.size)

if __name__ == "__main__":
    unittest.main()
