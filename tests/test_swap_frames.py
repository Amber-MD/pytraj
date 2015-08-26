from __future__ import print_function
import unittest; import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        from array import array
        list0 = [0, 1]
        list1 = [2, 3]
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]

        # a pair
        fa.swap(0, 1)
        aa_eq(fa[0].xyz, traj[1].xyz)
        aa_eq(fa[1].xyz, traj[0].xyz)

        # 2 pairs: from array
        fa = traj[:]
        fa.swap(array('i', list0), array('i', list1))
        aa_eq(fa[0].xyz, traj[2].xyz)
        aa_eq(fa[1].xyz, traj[3].xyz)

        # 2 pairs: from list
        fa = traj[:]
        fa.swap(list0, list1)
        aa_eq(fa[0].xyz, traj[2].xyz)
        aa_eq(fa[1].xyz, traj[3].xyz)


if __name__ == "__main__":
    unittest.main()
