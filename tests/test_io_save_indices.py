import unittest
from pytraj.base import *
import pytraj.io as mdio


class TestPyCpptrajIO(unittest.TestCase):

    def test_load_and_save_0(self):
        traj = mdio.iterload(filename="./data/Tc5b.x",
                             top="./data/Tc5b.top")[:10]
        indices = list(range(2, 3, 5)) + [3, 8, 9, 8]

        mdio.write_traj(filename="./output/test_io_saved_.x",
                        traj=traj,
                        top="./data/Tc5b.top",
                        frame_indices=indices,
                        overwrite=True)

        # check frames
        traj2 = mdio.iterload(filename="./output/test_io_saved_.x",
                              top="./data/Tc5b.top")
        assert traj2.n_frames == len(indices)


if __name__ == "__main__":
    unittest.main()
