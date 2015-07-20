from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.datafiles import load_cpptraj_output, tz2_ortho_trajin
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        cout = load_cpptraj_output(tz2_ortho_trajin + """
        rmsd first @CA perres range 2-7""")
        d = pt.rmsd_perres(traj, ref=0, mask='@CA', range='2-7')
        aa_eq(cout.values, d)

if __name__ == "__main__":
    unittest.main()
