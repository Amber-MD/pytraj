from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]

        pt.io.to_pickle(f0, './output/frame_pk.pk')
        f1 = pt.io.read_pickle('./output/frame_pk.pk')
        aa_eq(f0.xyz, f1.xyz)


if __name__ == "__main__":
    unittest.main()
