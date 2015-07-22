from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        cout = pt.datafiles.load_cpptraj_output("""
        parm ./data/Tc5b.top
        trajin ./data/md1_prod.Tc5b.x
        """, with_traj=True)

        traj2 = cout[0]
        aa_eq(traj.xyz, traj2.xyz)

if __name__ == "__main__":
    unittest.main()
