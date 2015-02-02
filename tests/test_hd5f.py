
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.trajs.Traj_HD5F import HD5F
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test
from pytraj.load_HD5F import load_hd5f
import numpy as np

# note : API will be changed. 
class Test(unittest.TestCase):
    def test_0(self):
        h5 = HD5F()
        # FIXME : *** Error in `python': double free or corruption
        farray = h5.load_toframearray("./data/ala2.h5")
        print (farray)
        print (farray[0])
        print (farray[0, 0])

    def test_1(self):
        traj = load_hd5f("./data/ala2.h5")
        print (traj)

if __name__ == "__main__":
    unittest.main()
