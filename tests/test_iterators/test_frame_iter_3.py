import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        for i, frame in enumerate(traj.iterframe(1, 6, 2)):
            pass

        assert i == 2

        for i, frame in enumerate(traj.iterframe(1, 5, 1)):
            pass

        assert i == 3

        for i, frame in enumerate(traj.iterframe(stop=8)):
            pass

        assert i == 7

        for i, frame in enumerate(traj.iterframe(start=7, stop=8)):
            pass

        assert i == 0

    def test_1(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        act = adict['distance']
        dslist = DatasetList()
        act.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act.setup(traj.top)

        for frame in traj.iterframe(stop=5):
            act.compute(frame)

        dslist = DatasetList()
        act2 = adict['distance']
        act2.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act2.setup(traj.top)
        act2.compute(traj.iterframe(stop=5))
        assert act2.n_frames == 5


if __name__ == "__main__":
    unittest.main()
