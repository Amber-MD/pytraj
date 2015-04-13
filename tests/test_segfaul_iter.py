from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj._shared_methods import _frame_iter_master
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        it = _frame_iter_master(traj)

        for idx, frame in enumerate(it):
            pass
            #print ("segmentation faul if uncommenting #traj")
            #traj[idx]

        fa = traj[:]
        for idx, frame in enumerate(fa):
            fa[idx]

if __name__ == "__main__":
    unittest.main()
