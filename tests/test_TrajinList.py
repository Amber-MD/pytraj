import os
import unittest
from pytraj.base import *
from pytraj.TrajinList import TrajinList

class TestTrajinList(unittest.TestCase):
    def test_0(self):
        topname = "./data/Tc5b.top"
        trajoutname = "./data/test.x"
        refilename = "./data/Tc5b.nat.crd"
        trajinname = "./data/md1_prod.Tc5b.x"
        toplist = TopologyList()
        toplist.add_parm(topname)
        toplist.info()
        
        top = toplist[0]
        
        #creat TrajinList instance
        trajininput= """
        reference Tc5b.nat.crd
        """
        
        argIn = ArgList(trajininput)
        trajlist = TrajinList()
        trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "1 3")
        trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "4 9 2")
        trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "5 7")
        trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "1 last")
        print(trajlist.max_frames)
        print (trajlist[0].n_frames)
        
        trajlist2 = TrajinList()
        trajlist2.add_traj("./data/Test_RemdTraj/rem.nc.000", 
                          "./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7",
                          "*")
        print (trajlist2.max_frames)

    def test_1(self):
        topname = "./data/Tc5b.top"
        top = Topology(topname)

        trajlist = TrajinList()
        trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "1 3")
        trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "4 9 2")
        trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "5 7")
        trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "1 last")
        print(trajlist.max_frames)
        print (trajlist[0].n_frames)

        trajlist.top = top.copy()
        print (trajlist.top)
        for traj in trajlist:
            pass
            #for frame in traj:
            #    print (frame)
        # FIXME: segmentation fault
        #print (traj[0].size)

if __name__ == "__main__":
    unittest.main()
