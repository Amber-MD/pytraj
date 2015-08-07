import unittest
#from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils.check_and_assert import is_word_in_class_name
from pytraj import calculate, adict
from pytraj.utils.Timer import Timer
from pytraj.plotting import simple_plot
from pytraj.utils import _import
import numpy as np


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d0 = calculate(adict['distance'], traj, ':2@CA :10@CA')[0]
        print(type(d0))
        print((d0.size))
        print(np.asarray(d0[:]))

        cppout = np.loadtxt(
            "./data/CAres2_CAres10.Tc5b.dat",
            skiprows=1).transpose()[1]
        print(cppout[:10])
        assert_almost_equal(d0[:], cppout[:d0.size], decimal=3)

        d1 = calculate('distance', traj, ':2@CA :10@CA')[0]
        assert_almost_equal(d1[:], cppout[:d1.size], decimal=3)

        with Timer() as t:
            trajlist = [traj, ]
        print("time to add traj to list = ", t.time_gap)

        with Timer() as t:
            d2 = calculate('distance', trajlist, ':2@CA :10@CA')[0]
        print("time to do actions = ", t.time_gap)

        with Timer() as t:
            for traj in trajlist:
                for frame in traj:
                    pass
        print("time to do looping all frames = ", t.time_gap)

        np.testing.assert_almost_equal(d2[:][:10], cppout[:10], decimal=3)

        has_plot, _plt = _import('matplotlib.pyplot')
        print(has_plot)
        if has_plot:
            print("pass")
            pass
            #plt = _plt.pyplot
            # plt.xlabel('snapshot #')
            #simple_plot(d2, 'ro')
            # plt.close()

    def test_2(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        from pytraj.common_actions import calc_distance
        d0 = calc_distance(traj, ":2@CA :10@CA")
        print(d0[:])

    def test_3(self):
        # load and calc_distance at the same time
        traj = mdio.load(*("./data/md1_prod.Tc5b.x", "./data/Tc5b.top"))
        from pytraj.common_actions import calc_distance
        d0 = calc_distance(traj, ":2@CA :10@CA", dtype='dataset')
        assert is_word_in_class_name(d0, 'Dataset')


if __name__ == "__main__":
    unittest.main()
