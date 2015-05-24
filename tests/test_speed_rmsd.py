from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
from pytraj.utils import Timer
import pytraj.common_actions as pyca

tip3p_dir = "./data/nogit/tip3p/"

class Test(unittest.TestCase):
    @test_if_path_exists(tip3p_dir)
    def test_0(self):
        fname = tip3p_dir  + "/md.trj"
        topname = tip3p_dir + "/tc5bwat.top"
        traj = mdio.iterload(fname, topname)
        ref = traj[10]
        print (traj)

        @Timer()
        def test_rmsd_pytraj_mode(traj, ref):
            traj.rmsd(ref=ref, mode='pytraj')

        @Timer()
        def test_rmsd_cpptraj_mode(traj, ref):
            traj.rmsd(ref=ref, mode='cpptraj')

        print ("test_rmsd_pytraj_mode")
        test_rmsd_pytraj_mode(traj, ref)
        print ("test_rmsd_cpptraj_mode")
        test_rmsd_cpptraj_mode(traj, ref)
         
if __name__ == "__main__":
    unittest.main()
