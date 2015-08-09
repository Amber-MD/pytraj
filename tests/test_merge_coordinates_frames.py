from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        # merge_coordinates
        import numpy as np

        # load 2 frames
        traj = pt.iterload(
            "./data/md1_prod.Tc5b.x", "./data/Tc5b.top",
            frame_slice=(0, 2))

        # test mutable traj
        coords = pt.tools.merge_coordinates(traj.to_mutable_trajectory())
        r0 = pt.misc.rmsd(coords, pt.get_coordinates(traj), True)
        assert r0 == 0.

        # test immutable traj
        coords = pt.tools.merge_coordinates(traj)
        r0 = pt.misc.rmsd(coords, pt.get_coordinates(traj), True)
        assert r0 == 0.

        # test tuple
        coords = pt.tools.merge_coordinates((frame for frame in traj))
        r0 = pt.misc.rmsd(coords, pt.get_coordinates(traj), True)
        assert r0 == 0.

        coords = pt.tools.merge_coordinates([f.copy() for f in traj])
        print(coords.shape)
        print(traj.xyz.shape)
        r0 = pt.misc.rmsd(coords, pt.get_coordinates(traj), True)
        print(r0)
        assert r0 == 0.

    def test_1(self):
        # merge_frames
        import numpy as np

        # load 2 frames
        traj = pt.iterload(
            "./data/md1_prod.Tc5b.x", "./data/Tc5b.top",
            frame_slice=(0, 2))

        # test mutable traj
        coords = pt.tools.merge_frames(traj.to_mutable_trajectory()).xyz
        r0 = pt.misc.rmsd(coords, pt.get_coordinates(traj), True)
        assert r0 == 0.

        #, True) test immutable traj
        assert np.any(pt.tools.merge_frames(traj).xyz.flatten() ==
                      pt.get_coordinates(traj).flatten())

        # tuple
        assert np.any(pt.tools.merge_frames((frame
                                             for frame in traj)).xyz.flatten()
                      == pt.get_coordinates(traj).flatten())

        # list
        assert np.any(pt.tools.merge_frames([frame
                                             for frame in traj]).xyz.flatten()
                      == pt.get_coordinates(traj).flatten())

        # frame_iter: all atoms
        assert np.any(pt.tools.merge_frames(traj()).xyz.flatten() ==
                      pt.get_coordinates(traj()).flatten())

        # frame_iter: CA atoms
        assert np.any(pt.tools.merge_frames(traj(mask='@CA')).xyz.flatten() ==
                      pt.get_coordinates(traj(mask='@CA')).flatten())


if __name__ == "__main__":
    unittest.main()
