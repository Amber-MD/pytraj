import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import calculate, adict
from pytraj.utils.Timer import Timer
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

        with Timer() as t:
            trajlist = [traj for _ in range(1000)]
        print ("time to add traj to list = ", t.time_gap)

        with Timer() as t:
            d2 = calculate('distance', ':2@CA :10@CA', trajlist)
        print ("time to do actions = ", t.time_gap)

        with Timer() as t:
            for traj in trajlist:
                for frame in traj:
                    pass
        print ("time to do looping all frames = ", t.time_gap)

        np.testing.assert_almost_equal(d2[:][:10], cppout[:10], decimal=3)
        print (d2.size)

if __name__ == "__main__":
    unittest.main()
