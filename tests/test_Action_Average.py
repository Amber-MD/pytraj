import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.common_actions import *
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

        # shorter
        from pytraj.common_actions import get_average_frame
        #frame2 = get_average_frame("", traj, traj.top)
        frame2 = get_average_frame("", traj)
        assert_almost_equal(frame2.coords, f_saved.coords)
        print (frame2[:2])
        print (f_saved[:2])
        print (traj[0, :2])
        print (frame2)
        print ("OK")

        frame3 = get_average_frame(traj=traj)
        assert_almost_equal(frame3.coords, f_saved.coords)

        # test list
        frame4 = get_average_frame(traj=[traj, traj[:3]], top=traj.top)
        print (frame4[:2])
        print (f_saved[:2])

        # test iter
        frame5 = get_average_frame(traj=traj(1, 7, 2), top=traj.top)
        f5_saved = mdio.load("./data/avg.Tc5b.frame_2_to_8_skip_2.pdb", traj.top)[0]
        assert_almost_equal(frame5.coords, f5_saved.coords)
        print (frame5[:2])
        print (f5_saved[:2])

        # test iter CA
        frame5 = get_average_frame("@CA", traj=traj(1, 7, 2), top=traj.top)
        print (frame5[:2])

        print ("frame5.n_atoms: for CA")
        print (frame5.n_atoms)

        # TODO: add cpptraj output here. For some reasons, I can not use 'average' with 
        # @CA mask in cpptraj. 
        #f5_saved = mdio.load("./data/avg.Tc5b.frame_2_to_8_skip_2.CA.pdb", traj.top)[0]
        #print (f5_saved[:2])
        #assert_almost_equal(frame5.coords, f5_saved.coords)

if __name__ == "__main__":
    unittest.main()
