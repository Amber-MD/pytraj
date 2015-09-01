from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")

        cout = pt.datafiles.load_cpptraj_output("""
        parm ./data/tz2.ortho.parm7
        trajin ./data/tz2.ortho.nc
        rms first nofit
        rms first mass
        """)
        #print(cout)

        aa_eq(pt.rmsd(traj, nofit=True), cout[0])
        aa_eq(pt.rmsd(traj, mass=True), cout[1])


if __name__ == "__main__":
    unittest.main()
