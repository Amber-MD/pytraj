import unittest
import numpy as np
from pytraj.base import *
from pytraj import adict
from pytraj import io as pt
from pytraj.common_actions import *
from pytraj.testing import aa_eq
from pytraj.common_actions import get_average_frame


class TestAverageFrame(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['average']
        dslist = DatasetList()
        act("* crdset s1", traj, dslist=dslist)
        d0 = dslist[0]

        frame = d0.get_frame()

        # make sure we DO reproducing cpptraj output
        f_saved = pt.iterload("./data/avg.Tc5b.pdb", traj.top)[0]
        aa_eq(frame.coords, f_saved.coords)

        # shorter
        from pytraj.common_actions import get_average_frame
        #frame2 = get_average_frame("", traj, traj.top)
        frame2 = get_average_frame(traj)
        aa_eq(frame2.coords, f_saved.coords)

        frame3 = get_average_frame(traj=traj)
        aa_eq(frame3.coords, f_saved.coords)

        # test list
        frame4 = get_average_frame(traj=[traj, traj[:3]], top=traj.top)

        # test iter
        frame5 = get_average_frame(traj=traj(1, 8, 2), top=traj.top)
        f5_saved = pt.iterload(
            "./data/avg.Tc5b.frame_2_to_8_skip_2.pdb", traj.top)[0]
        aa_eq(frame5.coords, f5_saved.coords)

        # test iter CA
        frame5 = get_average_frame(traj[[0, 3, 7]], '@CA', top=traj.top)

        # test frame_indices
        frame6 = get_average_frame(traj, mask='@CA', frame_indices=[0, 3, 7])
        aa_eq(frame5.xyz, frame6.xyz)

        xyz_0= pt.get_coordinates(traj(1, 8, 2))
        xyz_1= np.array([frame.xyz.copy() for frame in traj.iterframe(frame_indices=range(1, 8, 2))])
        aa_eq(xyz_0, xyz_1)

    def test_1(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        f0 = get_average_frame(fa, "@CA")
        f1 = get_average_frame(fa, "@CA")
        aa_eq(f0.xyz, f1.xyz)


if __name__ == "__main__":
    unittest.main()
