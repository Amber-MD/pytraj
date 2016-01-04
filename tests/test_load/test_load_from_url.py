from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestURL(unittest.TestCase):

    @unittest.skip('not reliable test')
    def test_url(self):
        try:
            import parmed as pmd
            url = "http://ambermd.org/tutorials/advanced/tutorial1/files/polyAT.pdb"
            traj0 = pt.load("./data/polyAT.pdb")
            traj = pt.load_parmed(url, as_traj=True)
            aa_eq(traj0.xyz, traj.xyz)

            traj1 = pt.load(url)
            aa_eq(traj1.xyz, traj.xyz)
        except ImportError:
            pass


if __name__ == "__main__":
    unittest.main()
