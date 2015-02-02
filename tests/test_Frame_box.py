import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
import numpy as np
from array import array

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        assert frame0.has_box() == False
        box = frame0.get_box()
        print(box)
        assert frame0.get_box().btype == 'nobox'
        print(box.bname)

        bview = frame0.boxview
        bview[3:] = np.asarray([109.471220634, 109.471220634, 109.471220634])
        print(frame0.has_box())
        print(frame0.get_box().btype)
        assert frame0.get_box().btype == 'truncoct'
        assert frame0.get_box().bname == 'Trunc. Oct.'
        print(frame0.get_box().help())

if __name__ == "__main__":
    unittest.main()
