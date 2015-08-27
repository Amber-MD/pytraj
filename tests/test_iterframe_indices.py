from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        t0 = traj[:]
        indices = range(3)

        d0 = pt.radgyr(traj._iterframe_indices(indices), top=traj.top)
        d1 = pt.radgyr(traj[indices])

        aa_eq(d0, d1)
        #print(d0)


if __name__ == "__main__":
    unittest.main()
