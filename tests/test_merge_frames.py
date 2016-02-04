from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        import numpy as np

        # load 2 frames
        traj = pt.iterload("./data/Tc5b.x",
                           "./data/Tc5b.top",
                           frame_slice=(0, 2))

        # test mutable traj
        big_frame = pt.tools.merge_frames(traj[:])
        assert pt.tools.rmsd(big_frame.xyz.flatten(),
                             pt.get_coordinates(traj).flatten()) < 1E-5

        # test immutable traj
        big_frame = pt.tools.merge_frames(traj)
        assert pt.tools.rmsd(big_frame.xyz.flatten(),
                             pt.get_coordinates(traj).flatten()) < 1E-5

        # test tuple
        big_frame = pt.tools.merge_frames((frame for frame in traj))
        assert pt.tools.rmsd(big_frame.xyz.flatten(),
                             pt.get_coordinates(traj).flatten()) < 1E-5

        # test list
        big_frame = pt.tools.merge_frames([frame.copy() for frame in traj])
        assert pt.tools.rmsd(big_frame.xyz.flatten(),
                             pt.get_coordinates(traj).flatten()) < 1E-5


if __name__ == "__main__":
    unittest.main()
