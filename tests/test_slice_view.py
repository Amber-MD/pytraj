from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.compat import zip

class Test(unittest.TestCase):
    def test_0(self):

        # create TrajectoryIter (readonly)
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # convert to Trajectory
        #fa = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]

        # frame view
        f0 =fa[0]
        aa_eq(f0.xyz, fa[0].xyz)
        # assign new values
        f0[0, 0] = 1000.
        assert f0[0, 0] == 1000.
        aa_eq(f0.xyz, fa[0].xyz)

        # sub-Trajectory view
        fa1 = fa[2:8:2]
        aa_eq(fa1.xyz, fa[2:8:2].xyz)
        fa1[0, 0] = 200.
        assert fa1[0, 0, 0] == 200.
        assert fa[2, 0, 0] == 200.
        aa_eq(fa1.xyz, fa[2:8:2].xyz)

        # sub-Trajectory view 2
        mylist = [1, 5, 8]
        fa2 = Trajectory()
        fa2.top = fa.top
        for index in mylist:
            fa2.append(fa[index], copy=False)
        aa_eq(fa2.xyz, fa[[1, 5, 8]].xyz)
        fa2[0, 0] = 500.
        assert fa2[0, 0, 0] == 500.
        assert fa[1, 0, 0] == 500.
        aa_eq(fa2.xyz, fa[[1, 5, 8]].xyz)

        # make sure we can make a copy
        facp = fa.copy()
        facp[0, 0, 0] = 501.
        assert (facp[0].rmsd(fa[0]) > 1.)
        assert facp[0, 0, 0] == 501.
        assert fa[0, 0, 0] != 501.

if __name__ == "__main__":
    unittest.main()
