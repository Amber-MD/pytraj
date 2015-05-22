from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.utils import Timer
from pytraj.compat import range # range = xrange in python2

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.testing import make_fake_traj

        # make a fake Trajectory with 1000 frames, 100K atoms
        traj = make_fake_traj(1000, 100000)
        print (traj)

        @Timer()
        def test_frame_iter_speed():
            for frame in traj:
                frame

        n_frames = traj.n_frames
        @Timer()
        def test_frame_indexing_speed():
            for i in range(n_frames):
                traj[i]

        n_frames = traj.n_frames
        @Timer()
        def test_frame_indexing_with_enumerate_speed():
            for i, frame in enumerate(traj):
                frame

        print ("test_frame_iter_speed")
        test_frame_iter_speed()
        print ("test_frame_indexing_speed")
        test_frame_indexing_speed()
        print ("test_frame_indexing_with_enumerate_speed")
        test_frame_indexing_with_enumerate_speed()


if __name__ == "__main__":
    unittest.main()
