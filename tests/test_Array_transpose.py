from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        c = pt.center_of_mass(traj, dtype='dataset')[0]
        ct = c.T
        aa_eq(ct.values, c.values.T)

        # view
        ct.values[0, 0] = 100.
        assert c.values.T[0, 0] == 100.
        print(c, ct)


if __name__ == "__main__":
    unittest.main()
