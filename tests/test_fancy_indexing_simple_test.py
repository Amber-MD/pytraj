import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
import numpy as np

class Test(unittest.TestCase):
    def test_0(self):
        # create FrameArray from Trajing_Single
        # TODO : add more assert
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(traj)
        arr0 = traj[:, :, :]
        print(arr0.shape)
        print(type(arr0))

if __name__ == "__main__":
    unittest.main()
