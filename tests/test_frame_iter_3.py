import unittest; import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        for i, frame in enumerate(traj.frame_iter(1, 6, 2)):
            pass

        assert i == 2

        for i, frame in enumerate(traj.frame_iter(1, 5, 1)):
            pass

        assert i == 3

        for i, frame in enumerate(traj.frame_iter(stop=8)):
            pass

        assert i == 7

        for i, frame in enumerate(traj.frame_iter(start=7, stop=8)):
            #print(frame)

        assert i == 0

    def test_1(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['distance']
        dslist = DataSetList()
        #d0 = act(":2@CA :10@CA", (traj.frame_iter(stop=5),), traj.top, quick_get=True)
        act.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act.process(traj.top)
        #print(list(traj.frame_iter()))
        #print(list(traj.chunk_iter(4)))

        for frame in traj.frame_iter(stop=5):
            act.do_action(frame)
        #print(dslist.size)
        #print(dslist[0][:])

        dslist = DataSetList()
        act2 = adict['distance']
        act2.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act2.process(traj.top)
        act2.do_action(traj.frame_iter(stop=5))
        #print(act2.n_frames)
        #print(dslist.size)
        #print(dslist[0].size)
        assert act2.n_frames == 5
        #print(dslist[0][:])


if __name__ == "__main__":
    unittest.main()
