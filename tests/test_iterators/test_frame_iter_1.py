import unittest
import pytraj as pt
from pytraj.utils.check_and_assert import assert_almost_equal


class TestIterator(unittest.TestCase):

    def test_frame_iterator(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        frame = traj[0]

        for arr0 in frame:
            arr0[0] = 0.0

        # make sure that we update frame too
        for i in range(frame.n_atoms):
            assert (frame[i, 0]) == 0.0

        for f0 in pt.iterframe(frame):
            pass

        for f0 in pt.iterframe(frame):
            pass

        for frame in traj.iterframe():
            pass

        for frame in traj.iterframe(stop=traj.n_frames - 5):
            pass

        for frame in traj.iterframe(start=5, stop=traj.n_frames - 2):
            pass

        for frame in traj.iterframe(start=1, stop=traj.n_frames - 2, step=2):
            pass


if __name__ == "__main__":
    unittest.main()
