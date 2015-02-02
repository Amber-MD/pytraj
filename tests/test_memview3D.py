import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
import numpy as np

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
        arr0 = traj[:, :, :]
        print(arr0.shape)
        arr0[0, 0, 0] = 105.
        print(traj)
        print(traj[0, 0, 0])
        # TODO : return (n_frames, n_atoms, 3) array as a view of traj

if __name__ == "__main__":
    unittest.main()
