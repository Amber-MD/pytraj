from __future__ import print_function
import unittest

from pytraj import io as pt
from pytraj.testing import aa_eq

from utils import fn


class Test(unittest.TestCase):
    def test_0(self):

        # create TrajectoryIter (readonly)
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        # convert to Trajectory
        fa = traj[:]

        # frame view
        f0 = fa[0]
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

        # make sure we can make a copy
        facp = fa.copy()
        facp.xyz[0, 0, 0] = 501.
        assert facp.xyz[0, 0, 0] == 501.
        assert fa.xyz[0, 0, 0] != 501.


if __name__ == "__main__":
    unittest.main()
