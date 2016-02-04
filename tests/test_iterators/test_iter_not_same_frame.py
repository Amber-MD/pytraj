from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        # TrajectoryIterator
        aa_eq(
            pt.molsurf([f.copy() for f in traj],
                       top=traj.top),
            pt.molsurf(traj))

        # frame_iter
        aa_eq(
            pt.molsurf([f.copy() for f in traj()],
                       top=traj.top),
            pt.molsurf(traj))

        aa_eq(
            pt.molsurf([f.copy() for f in traj.iterchunk(3)],
                       top=traj.top),
            pt.molsurf(traj))


if __name__ == "__main__":
    unittest.main()
