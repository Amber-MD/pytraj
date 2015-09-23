import unittest
import numpy as np
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj.common_actions import *
from pytraj.testing import aa_eq
from pytraj.common_actions import mean_structure
from pytraj import Trajectory


class TestAverageFrame(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        # make sure we DO reproducing cpptraj output
        f_saved = pt.iterload("./data/avg.Tc5b.pdb", traj.top)[0]

        # shorter
        from pytraj.common_actions import mean_structure
        #frame2 = mean_structure("", traj, traj.top)
        frame2 = mean_structure(traj)
        aa_eq(frame2.coords, f_saved.coords, decimal=3)

        frame3 = mean_structure(traj=traj)
        aa_eq(frame3.coords, f_saved.coords, decimal=3)

        # test list
        frame4 = mean_structure(traj=[traj, traj[:3]], top=traj.top)

        # test iter
        frame5 = mean_structure(traj=traj(1, 8, 2), top=traj.top)
        f5_saved = pt.iterload(
            "./data/avg.Tc5b.frame_2_to_8_skip_2.pdb", traj.top)[0]
        aa_eq(frame5.coords, f5_saved.coords, decimal=3)

        # test iter CA
        frame5 = mean_structure(traj[[0, 3, 7]], '@CA', top=traj.top)

        # test frame_indices
        frame6 = mean_structure(traj, mask='@CA', frame_indices=[0, 3, 7])
        aa_eq(frame5.xyz, frame6.xyz, decimal=3)

        xyz_0= pt.get_coordinates(traj(1, 8, 2))
        xyz_1= np.array([frame.xyz.copy() for frame in traj.iterframe(frame_indices=range(1, 8, 2))])
        aa_eq(xyz_0, xyz_1, decimal=3)

        # test as traj
        out_traj = mean_structure(traj, mask='@CA', frame_indices=[0, 3, 7], restype='traj')
        assert isinstance(out_traj, Trajectory), 'must be Trajectory'
        aa_eq(out_traj.xyz, frame6.xyz, decimal=3)


    def test_1(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        f0 = mean_structure(fa, "@CA")
        f1 = mean_structure(fa, "@CA")


if __name__ == "__main__":
    unittest.main()
