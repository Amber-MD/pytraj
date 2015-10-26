from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca


class TestAutoImageAndRotateDihedral(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        farray = traj[:]

        t0trajectory = pt.trajectory.Trajectory(traj)
        aa_eq(farray.unitcells, t0trajectory.unitcells)

        # autoimage
        farray.autoimage()
        t0trajectory.autoimage()
        aa_eq(farray.xyz, t0trajectory.xyz)

        # rotate_dihedral
        pt.rotate_dihedral(t0trajectory, '3:phi:120')
        pt.rotate_dihedral(farray, '3:phi:120')
        aa_eq(farray.xyz, t0trajectory.xyz)

        aa_eq(pt.calc_phi(t0trajectory, '3').values, [120
              for _ in range(t0trajectory.n_frames)])


class TestNoName(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        trajectory_traj = traj[:]

        # test xyz
        aa_eq(trajectory_traj.xyz, traj.xyz)

        # test object lifetime
        aa_eq(trajectory_traj[0].xyz, trajectory_traj.xyz[0])

        # test Box
        assert (trajectory_traj.top.has_box() == True)
        boxes = traj.unitcells
        for i, frame in enumerate(trajectory_traj):
            assert (frame.has_box() == True)
            f_blist = frame.box.tolist()
            aa_eq(f_blist, boxes[i].tolist())

        # test autoimage
        # make Trajectory from TrajectoryIterator
        fa = traj[:]
        fa.autoimage()
        saved_traj = pt.iterload(
            "./data/tz2.autoimage.nc", "./data/tz2.ortho.parm7")

        # make sure to reproduce cpptraj's output too
        aa_eq(saved_traj.xyz, fa.xyz)

    def testFromIterable(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        aa_eq(pt.trajectory.Trajectory.from_iterable(traj).xyz, traj.xyz)


if __name__ == "__main__":
    unittest.main()
