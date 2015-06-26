from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        from glob import glob
        fname = "./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7"
        flist = sorted(glob("./data/Test_RemdTraj/rem.nc.*")) 
        print (flist)

        traj = pt.iterload(flist, fname)

        # 4 chunks
        trajlist = list(traj.split_iterators(4))
        aa_eq(trajlist[0].xyz, traj[:10].xyz)
        aa_eq(trajlist[1].xyz, traj[10:20].xyz)
        aa_eq(trajlist[2].xyz, traj[20:30].xyz)
        aa_eq(trajlist[3].xyz, traj[30:].xyz)

        # 2 chunks
        trajlist = list(traj.split_iterators(2))
        aa_eq(trajlist[0].xyz, traj[:20].xyz)
        aa_eq(trajlist[1].xyz, traj[20:].xyz)

if __name__ == "__main__":
    unittest.main()
