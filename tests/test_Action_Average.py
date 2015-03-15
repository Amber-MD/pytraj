import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['average']
        dslist = DataSetList()
        act("* crdset s1", traj, dslist=dslist)
        act.print_output()
        print (dslist.size)
        d0 = dslist[0]
        print (d0.dtype)
        print (d0.name)

        frame = d0.get_frame()
        print (frame)

        # make sure we DO reproducing cpptraj output
        f_saved = mdio.load("./data/avg.Tc5b.pdb", traj.top)[0]
        assert_almost_equal(frame.coords, f_saved.coords)
        print ("OK")

if __name__ == "__main__":
    unittest.main()
