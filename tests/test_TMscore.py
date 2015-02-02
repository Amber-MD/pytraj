import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.tm_score import tm_score

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # we could generate new traj.top from pdb file too
        #traj.top = Topology("./test_0_after.pdb")
        tm_0 = tm_score(traj[0], traj[1])
        #print tm_0
        #print tm_score(traj[0], traj[0]) # == 1.
        #print tm_score(traj[0], traj[1], mask="!@CA", top=traj.top)
        for idx, frame in enumerate(traj[490:501] + traj[-1:-10:-1]):
            score = tm_score(frame, traj[0], mask1="!:2-10@CA", top1=traj.top,
                                             mask2="!:2-10@CA", top2=traj.top)
            print("%s : %s" % (idx, score))

if __name__ == "__main__":
    unittest.main()
