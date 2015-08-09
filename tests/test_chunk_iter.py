import unittest
from pytraj import *
from pytraj.base import *
from pytraj.common_actions import calc_distance
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal, is_generator
from pytraj.decorators import no_test
import numpy as np


class Test(unittest.TestCase):
    #@no_test

    def test_0(self):
        dslist = DataSetList()
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        calc_distance = adict['distance']

        act = calc_distance
        act.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act.process(traj.top)
        act.do_action((traj.chunk_iter(2), ))
        assert dslist.size == 1
        print(dslist[0].size)
        assert dslist[0].size == traj.n_frames
        cppout = np.loadtxt(
            "./data/CAres2_CAres10.Tc5b.dat",
            skiprows=1).transpose()[1]
        assert_almost_equal(dslist[0][:], cppout)
        act.do_action((traj.chunk_iter(chunksize=4, stop=8), ))
        print(act.n_frames)
        print(dslist[0].size)

    def test_1(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        i = 0
        for farray in traj.chunk_iter(chunksize=4, stop=8):
            print(farray)
            i += 1
            print(i)
        assert i == 2
        print("n-chunk = ", i)

        i = 0
        for farray in traj.chunk_iter(chunksize=4):
            print(farray)
            i += 1
            print(i)
        assert i == 3
        print("n-chunk = ", i)

        i = 0
        for farray in traj.chunk_iter(chunksize=2):
            print(farray)
            i += 1
            print(i)
        assert i == 5
        print("n-chunk = ", i)

        print("for farray in traj.chunk_iter(start=3, chunksize=4, stop=8)")
        i = 0
        for farray in traj.chunk_iter(start=3, chunksize=4, stop=8):
            print(farray)
            i += 1
            print(i)
        assert i == 2
        print("n-chunk = ", i)

        # action on chunk_iter
        import pytraj.common_actions as pyca
        pyca.calc_distance(
            [traj.chunk_iter(), traj.chunk_iter(), traj[0]], '@CA @CB',
            top=traj.top)


if __name__ == "__main__":
    unittest.main()
