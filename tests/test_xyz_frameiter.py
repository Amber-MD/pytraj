from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        aa_eq(traj(stop=3).xyz, traj[:3].xyz)
        aa_eq(traj(stop=3, mask='@CA').xyz, traj[:3, '@CA'].xyz)

if __name__ == "__main__":
    unittest.main()
