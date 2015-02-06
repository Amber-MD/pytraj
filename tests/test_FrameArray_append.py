import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = FrameArray()
        for i in range(1000):
            for frame in traj:
                farray.append(frame, copy=False)

        print (farray)

        farray.join(farray)
        print (farray)

if __name__ == "__main__":
    unittest.main()
