import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import calculate, adict
import numpy as np

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d0 = calculate(adict['distance'], ':2@CA :10@CA', traj)
        print (d0.size)

        cppout = np.loadtxt("./data/CAres2_CAres10.Tc5b.dat", skiprows=1).transpose()[1]
        print(cppout[:10])
        np.testing.assert_almost_equal(d0[:], cppout[:d0.size], decimal=3)

        d1 = calculate('distance', ':2@CA :10@CA', traj)
        np.testing.assert_almost_equal(d1[:], cppout[:d1.size], decimal=3)

if __name__ == "__main__":
    unittest.main()
