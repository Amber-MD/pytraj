import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_1(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        act2 = adict['distance']
        act2.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act2.process(traj.top)
        act2.do_action(traj.chunk_iter())
        assert act2.n_frames == 10
        print(dslist[0][:])

if __name__ == "__main__":
    unittest.main()
