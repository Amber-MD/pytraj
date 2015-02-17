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
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        calc_distance = adict['distance']

        act = calc_distance
        act.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act.process(traj.top)
        act.do_action(0, (traj.chunk_iter(2),))
        assert dslist.size == 1
        print (dslist[0].size)
        assert dslist[0].size == traj.n_frames
        cppout = np.loadtxt("./data/CAres2_CAres10.Tc5b.dat", skiprows=1).transpose()[1]
        assert_almost_equal(dslist[0][:], cppout)
        act.do_action(0, (traj.chunk_iter(chunk=4, stop=8),))
        print (dslist[0].size)

        # TODO : fail
        #assert dslist[0].size == traj.n_frames + 8

    def test_1(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        i = 0
        for farray in traj.chunk_iter(chunk=4, stop=8):
            print (farray)
            i += 1
            print (i)
        assert i == 2
        print ("n-chunk = ", i)

        i = 0
        for farray in traj.chunk_iter(chunk=4):
            print (farray)
            i += 1
            print (i)
        assert i == 3
        print ("n-chunk = ", i)

        i = 0
        for farray in traj.chunk_iter(chunk=2):
            print (farray)
            i += 1
            print (i)
        assert i == 5
        print ("n-chunk = ", i)

        print ("for farray in traj.chunk_iter(start=3, chunk=4, stop=8)")
        i = 0
        for farray in traj.chunk_iter(start=3, chunk=4, stop=8):
            print (farray)
            i += 1
            print (i)
        assert i == 2
        print ("n-chunk = ", i)

if __name__ == "__main__":
    unittest.main()
