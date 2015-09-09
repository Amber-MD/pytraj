import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.common_actions import *
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.testing import aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['average']
        dslist = DatasetList()
        act("* crdset s1", traj, dslist=dslist)
        d0 = dslist[0]

        frame = d0.get_frame()

        # make sure we DO reproducing cpptraj output
        f_saved = mdio.iterload("./data/avg.Tc5b.pdb", traj.top)[0]
        assert_almost_equal(frame.coords, f_saved.coords)

        # shorter
        from pytraj.common_actions import get_average_frame
        #frame2 = get_average_frame("", traj, traj.top)
        frame2 = get_average_frame(traj)
        assert_almost_equal(frame2.coords, f_saved.coords)

        frame3 = get_average_frame(traj=traj)
        assert_almost_equal(frame3.coords, f_saved.coords)

        # test list
        frame4 = get_average_frame(traj=[traj, traj[:3]], top=traj.top)

        # test iter
        frame5 = get_average_frame(traj=traj(1, 8, 2), top=traj.top)
        f5_saved = mdio.iterload(
            "./data/avg.Tc5b.frame_2_to_8_skip_2.pdb", traj.top)[0]
        assert_almost_equal(frame5.coords, f5_saved.coords)

        # test iter CA
        frame5 = get_average_frame(traj(1, 8, 2), '@CA', top=traj.top)

        # TODO: add cpptraj output here. For some reasons, I can not use 'average' with
        # @CA mask in cpptraj.
        #f5_saved = mdio.iterload("./data/avg.Tc5b.frame_2_to_8_skip_2.CA.pdb", traj.top)[0]
        #assert_almost_equal(frame5.coords, f5_saved.coords)

    def test_1(self):
        from pytraj.utils import Timer
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        from pytraj.common_actions import get_average_frame
        f0 = get_average_frame(fa, "@CA")
        f1 = get_average_frame(fa, "@CA")
        aa_eq(f0.xyz, f1.xyz)


def average_cpptraj(fa):
    get_average_frame(fa)


def average_pytraj(fa):
    fa.average()


if __name__ == "__main__":
    unittest.main()
