import unittest
from pytraj import *
from pytraj.base import *
from pytraj.common_actions import calc_distance
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal, is_generator

class Test(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        calc_distance = adict['distance']

        #calc_distance(":2@CA :10@CA", traj.chunk_iter(2), traj.top, dslist=dslist)
        act = calc_distance
        act.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act.process(traj.top)
        act.do_action(0, traj.chunk_iter(2))
        #for i, farray in enumerate(traj.chunk_iter(2)):
        #    act.do_action(0, farray)
        print (dslist.size)
        print (dslist[0].size)

if __name__ == "__main__":
    unittest.main()
