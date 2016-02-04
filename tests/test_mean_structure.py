#!/usr/bin/env python
import unittest
import numpy as np
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj.all_actions import *
from pytraj.testing import aa_eq
from pytraj import Trajectory
from pytraj import mean_structure


class TestAverageFrame(unittest.TestCase):

    def test_comprehensive(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        # make sure we DO reproducing cpptraj output
        f_saved = pt.iterload("./data/avg.Tc5b.pdb", traj.top)[0]

        # shorter
        frame2 = mean_structure(traj)
        aa_eq(frame2.xyz, f_saved.xyz, decimal=3)

        frame3 = mean_structure(traj=traj)
        aa_eq(frame3.xyz, f_saved.xyz, decimal=3)

        # test list
        frame4 = mean_structure(traj=[traj, traj[:3]], top=traj.top)

        # test iter
        frame5 = mean_structure(traj=traj(1, 8, 2), top=traj.top)
        f5_saved = pt.iterload("./data/avg.Tc5b.frame_2_to_8_skip_2.pdb",
                               traj.top)[0]
        aa_eq(frame5.xyz, f5_saved.xyz, decimal=3)

        # test iter CA
        frame5 = mean_structure(traj[[0, 3, 7]], '@CA', top=traj.top)

        # use atom_indices
        ca_indices = pt.select('@CA', traj.top)
        frame5_1 = mean_structure(traj[[0, 3, 7]], ca_indices, top=traj.top)

        # test frame_indices
        frame6 = mean_structure(traj, mask='@CA', frame_indices=[0, 3, 7])
        aa_eq(frame5.xyz, frame6.xyz, decimal=3)
        aa_eq(frame5_1.xyz, frame6.xyz, decimal=3)

        xyz_0 = pt.get_coordinates(traj(1, 8, 2))
        xyz_1 = np.array([frame.xyz.copy(
        ) for frame in traj.iterframe(frame_indices=range(1, 8, 2))])
        aa_eq(xyz_0, xyz_1, decimal=3)

        # test as traj
        out_traj = mean_structure(traj,
                                  mask='@CA',
                                  frame_indices=[0, 3, 7],
                                  dtype='traj')
        assert isinstance(out_traj, Trajectory), 'must be Trajectory'
        aa_eq(out_traj.xyz, frame6.xyz, decimal=3)

        # raise if not trajectory, traj or frame
        self.assertRaises(ValueError, lambda: pt.mean_structure(traj, dtype='trajxyz'))

    def test_autoimage(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        t0 = traj[:]

        t0.autoimage()
        avg_0 = pt.mean_structure(t0, '@CA')
        avg_1 = pt.mean_structure(traj(autoimage=True), '@CA')
        aa_eq(avg_0.xyz, avg_1.xyz)

    def test_autoimage_with_rmsfit(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        t0 = traj[:]

        pt.autoimage(t0).superpose()
        avg_0 = pt.mean_structure(t0, '@CA')
        avg_1 = pt.mean_structure(traj(autoimage=True, rmsfit=0), '@CA')
        aa_eq(avg_0.xyz, avg_1.xyz)

        # 3rd frame
        # assign traj again
        t0 = traj[:]
        pt.autoimage(t0).superpose(ref=3)
        avg_0 = pt.mean_structure(t0, '@CA')
        avg_1 = pt.mean_structure(traj(autoimage=True, rmsfit=3), '@CA')
        avg_2 = pt.mean_structure(traj, autoimage=True, rmsfit=3, mask='@CA')
        aa_eq(avg_0.xyz, avg_1.xyz)
        aa_eq(avg_0.xyz, avg_2.xyz)

        # 3rd frame, frame_indices
        # assign traj again
        frame_indices = [0, 8, 5]
        t0 = traj[frame_indices]
        t1 = traj[frame_indices]

        t0.autoimage().superpose(ref=-1)
        avg_0 = pt.mean_structure(t0, '@CA')

        # use ref=5 which correspond to original index
        # try with pytraj.TrajectoryIterator
        avg_1 = pt.mean_structure(traj,
                                  autoimage=True,
                                  rmsfit=5,
                                  mask='@CA',
                                  frame_indices=frame_indices)
        # try with pytraj.Trajectory
        avg_2 = pt.mean_structure(t1, autoimage=True, rmsfit=-1, mask='@CA')
        avg_3 = pt.mean_structure(traj[:],
                                  autoimage=True,
                                  rmsfit=5,
                                  mask='@CA',
                                  frame_indices=frame_indices)

        aa_eq(avg_0.xyz, avg_1.xyz)
        aa_eq(avg_0.xyz, avg_2.xyz)
        aa_eq(avg_0.xyz, avg_3.xyz)


if __name__ == "__main__":
    unittest.main()
