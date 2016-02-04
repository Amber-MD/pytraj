from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestIteraframeIndices(unittest.TestCase):

    def test_iterframe_indices(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        t0 = traj[:]
        indices = range(3)

        d0 = pt.radgyr(traj._iterframe_indices(indices), top=traj.top)
        d1 = pt.radgyr(traj[indices])
        d2 = pt.radgyr(traj, frame_indices=indices)

        aa_eq(d0, d1)
        aa_eq(d0, d2)

        # raise if out of bound
        # only care about TrajectoryCpptraj since we would get segmentation fault
        # if index is larger than max n_frame
        def iter_(traj=traj):
            for _ in traj._iterframe_indices([traj.n_frames, ]):
                print(_.xyz)

        self.assertRaises(AssertionError, lambda: iter_(traj))


if __name__ == "__main__":
    unittest.main()
