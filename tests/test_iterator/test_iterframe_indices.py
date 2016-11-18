from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.utils import aa_eq
import pytest


class TestIteraframeIndices(unittest.TestCase):

    def test_iterframe_indices(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        t0 = traj[:]
        indices = range(3)

        d1 = pt.radgyr(traj[indices])
        d2 = pt.radgyr(traj, frame_indices=indices)

        aa_eq(d2, d1)

        # raise if out of bound
        # only care about TrajectoryCpptraj since we would get segmentation fault
        # if index is larger than max n_frame
        with pytest.raises(AssertionError):
            for _ in traj._iterframe_indices([traj.n_frames, ]):
                print(_.xyz)

if __name__ == "__main__":
    unittest.main()
