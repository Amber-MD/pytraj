from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.parallel.base import WrapBareIterator


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        # TrajectoryIterator
        myiter = WrapBareIterator([f.copy() for f in traj], top=traj.top)
        aa_eq(pt.molsurf(myiter), pt.molsurf(traj))

        # frame_iter
        myiter = WrapBareIterator([f.copy() for f in traj()], top=traj.top)
        aa_eq(pt.molsurf(myiter), pt.molsurf(traj))

        myiter = WrapBareIterator(
            [f.copy() for f in traj.iterchunk(3)], top=traj.top)
        aa_eq(pt.molsurf(myiter), pt.molsurf(traj))


if __name__ == "__main__":
    unittest.main()
