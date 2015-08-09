from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        url = "http://ambermd.org/tutorials/advanced/tutorial1/files/polyAT.pdb"
        traj0 = pt.load_pdb("./data/polyAT.pdb")
        traj = pt.load_parmed(url, as_traj=True)
        aa_eq(traj0.xyz, traj.xyz)

        traj1 = pt.load(url)
        aa_eq(traj1.xyz, traj.xyz)


if __name__ == "__main__":
    unittest.main()
