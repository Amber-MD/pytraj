import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.common_actions import *
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.testing import aa_eq
from pytraj.testing import make_fake_traj

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
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
        f_saved = mdio.iterload("./data/avg.Tc5b.pdb", traj.top)[0]
        assert_almost_equal(frame.coords, f_saved.coords)

        # shorter
        from pytraj.common_actions import get_average_frame
        #frame2 = get_average_frame("", traj, traj.top)
        frame2 = get_average_frame(traj)
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
        f5_saved = mdio.iterload("./data/avg.Tc5b.frame_2_to_8_skip_2.pdb", traj.top)[0]
        assert_almost_equal(frame5.coords, f5_saved.coords)
        print (frame5[:2])
        print (f5_saved[:2])

        # test iter CA
        frame5 = get_average_frame(traj(1, 7, 2), '@CA', top=traj.top)
        print (frame5[:2])

        print ("frame5.n_atoms: for CA")
        print (frame5.n_atoms)

        # TODO: add cpptraj output here. For some reasons, I can not use 'average' with 
        # @CA mask in cpptraj. 
        #f5_saved = mdio.iterload("./data/avg.Tc5b.frame_2_to_8_skip_2.CA.pdb", traj.top)[0]
        #print (f5_saved[:2])
        #assert_almost_equal(frame5.coords, f5_saved.coords)

    def test_1(self):
        from pytraj.utils import Timer
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        from pytraj.common_actions import get_average_frame
        f0 = get_average_frame(fa, "@CA")
        f1 = fa.average("@CA")
        print (f0.rmsd_nofit(f1))
        aa_eq(f0.xyz, f1.xyz)

        f0 = get_average_frame(fa, "!@H=")
        f1 = fa.average("!@H=")
        print (f0.rmsd_nofit(f1))
        aa_eq(f0.xyz, f1.xyz)

        f0 = get_average_frame(fa)
        f1 = fa.average()
        print (f0.rmsd_nofit(f1))
        aa_eq(f0.xyz, f1.xyz)

        @Timer()
        def average_pytraj(fa):
            fa.average()

        @Timer()
        def average_cpptraj(fa):
            get_average_frame(fa)

        print (fa)
        print ("average_pytraj")
        average_pytraj(fa)
        print ("average_cpptraj")
        average_cpptraj(fa)

        from pytraj.testing import make_fake_traj
        fa = make_fake_traj(100, 10000)
        print (fa)
        print ("average_pytraj")
        average_pytraj(fa)
        print ("average_cpptraj")
        average_cpptraj(fa)

def average_cpptraj(fa):
    get_average_frame(fa)

def average_pytraj(fa):
    fa.average()

if __name__ == "__main__":
    unittest.main()
    #from memory_profiler import memory_usage
    #from numpy import max
    #fa = make_fake_traj(10000, 10000)
    #m_pytraj = max(memory_usage((average_pytraj, (fa,))))
    #m_cpptraj = max(memory_usage((average_cpptraj, (fa,))))
    #print (m_cpptraj, m_pytraj)
