import unittest
from pytraj import io as mdio
import numpy as np


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        assert frame0.has_box() == False
        frame0.box
        assert frame0.box.type == 'nobox'

        bview = frame0._boxview
        bview[3:] = np.asarray([109.471220634, 109.471220634, 109.471220634])
        bview[:3] = np.array([100, 100, 100], dtype='f8')
        assert frame0.box.type == 'truncoct'


if __name__ == "__main__":
    unittest.main()
