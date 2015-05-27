import unittest
from pytraj.compat import set
from glob import glob
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
       flist = glob("./data/Test_RemdTraj/rem.nc.*") 
       trajlist = []
       for fh in flist:
           topfile = "./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7"
           trajlist.append(mdio.iterload(fh, topfile))

       Tset = set([])
       f4922 = Trajectory()
       f4922.resize(trajlist[0].n_frames)
       print(f4922.n_frames)
       f4922.top = trajlist[0].top.copy()

       for traj in trajlist:
           for idx, frame in enumerate(traj):
               if frame.temperature == 492.2:
                   f4922[idx] = frame

       print(f4922.temperatures)
       print(f4922[0, 0, :])

       # make sure we reproduce cpptraj output
       cpptraj = mdio.iterload("./data/Test_RemdTraj/temp0.crd.492.20", topfile)

       for idx, framepy in enumerate(f4922):
           assert_almost_equal(framepy.coords, cpptraj[idx].coords)

       print("YES")

if __name__ == "__main__":
    unittest.main()
