import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
import numpy as np
from array import array


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        assert frame0.has_box() == False
        box = frame0.box
        print(box)
        assert frame0.box.type == 'nobox'
        print(box.name)

        bview = frame0.boxview
        bview[3:] = np.asarray([109.471220634, 109.471220634, 109.471220634])
        print(frame0.has_box())
        assert frame0.box.type == 'truncoct'
        assert frame0.box.name == 'Trunc. Oct.'


if __name__ == "__main__":
    unittest.main()
