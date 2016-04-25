#!/usr/bin/env python
import unittest
import pytraj as pt
from pytraj import adict
from pytraj.testing import aa_eq
from pytraj.utils import tempfolder


class TestSuperposeTrajectory(unittest.TestCase):
    """superpose mutable Trajectory
    """

    def test_frame_fit(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        f1 = traj[1]

        arr0 = list(f0[0])
        arr1 = list(f1[0])

        f0.rmsd(f1)
        aa_eq(arr0, f0[0])
        aa_eq(arr1, f1[0])

        f1.rmsfit(ref=f0)

        # expect reference `f0` xyz are not changed
        aa_eq(arr0, f0[0])

        trajsaved = pt.iterload("./data/fit_to_1stframe.Tc5b.x",
                                "./data/Tc5b.top")
        f1saved = trajsaved[1]

        # make sure we reproduce cpptraj output
        aa_eq(f1.xyz, f1saved.xyz, decimal=3)

        farray = traj[:]
        farray.rmsfit(ref=traj[0])
        aa_eq(farray[1].xyz, f1saved.xyz, decimal=3)

    def test_0(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        f0 = traj[0]
        f0saved = f0.copy()
        f1 = traj[1]

        rmsd_0 = f0.rmsd(f1)
        rmsd_0_nofit = f0.rmsd_nofit(f1)
        assert rmsd_0 != rmsd_0_nofit

        # do fitting
        f1.rmsfit(ref=f0)
        rmsd_1 = f1.rmsd(f0)
        rmsd_1_nofit = f1.rmsd_nofit(f0)

        # make sure that rmsd_nofit after do fitting is equal to rmsd (with
        # fit)
        assert rmsd_1 - rmsd_1_nofit < 1E-3

        farray.rmsfit(ref=f0)
        assert rmsd_1 - farray[1].rmsd_nofit(f0) < 1E-3

    def test_1(self):

        # load frames to immutable traj
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        trajsaved = pt.iterload("./data/fit_to_1stframe.Tc5b.x",
                                "./data/Tc5b.top")

        for _f1 in trajsaved:
            pass

        f0saved = traj[0].copy()
        first = traj[0].copy()

        # make mutable traj
        farray = traj[:]

        aa_eq(farray[0].xyz, first.xyz)
        farray.rmsfit(ref=first, mask="*", mass=False)
        farray2 = traj[:]
        farray2.superpose(ref=first, mask="*", mass=False)

        for i, _f0 in enumerate(farray):
            _f1 = trajsaved[i]
            aa_eq(_f0.xyz, _f1.xyz, decimal=3)

        for i, _f0 in enumerate(farray2):
            _f1 = trajsaved[i]
            aa_eq(_f0.xyz, _f1.xyz, decimal=3)

    def test_frame_indices(self):
        # load frames to immutable traj
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        # convert to numpy traj

        frame_indices = [0, 3, 5]

        t00 = traj[:]
        t01 = traj[:]
        t10 = traj[frame_indices]
        t11 = traj[frame_indices]
        aa_eq(t00[frame_indices].xyz, t11.xyz)

        ref = traj[-1]
        t00.superpose(ref=ref, frame_indices=frame_indices)

        ref = traj[-1]
        pt.superpose(t01, ref=ref, frame_indices=frame_indices)

        ref = traj[-1]
        t10.superpose(ref=ref)

        ref = traj[-1]
        pt.superpose(t11, ref=ref, frame_indices=frame_indices)

        aa_eq(t00.xyz, t01.xyz)
        aa_eq(t10.xyz, t11.xyz)
        aa_eq(t00[frame_indices].xyz, t10.xyz)

    def testsuperpose_vs_rmsd(self):
        # load frames to immutable traj
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        t0 = traj[:]
        t1 = traj[:]
        pt.rmsd(t0, ref=traj[0], mask='@CA')
        pt.superpose(t1, ref=traj[0], mask='@CA')
        aa_eq(t0.xyz, t1.xyz)

class TestSuperposeTrajectoryIterator(unittest.TestCase):
    """test superpose TrajectoryIterator
    """

    def test_superpose_trajectory_iterator(self):
        traj_immut = pt.iterload("data/Tc5b.x", "data/Tc5b.top")
        traj_immut2 = pt.iterload("data/Tc5b.x", "data/Tc5b.top")
        traj_mut = pt.load("data/Tc5b.x", "data/Tc5b.top")

        ref = pt.iterload("data/Tc5b.crd", "data/Tc5b.top")[0]
        traj_mut.superpose(ref=ref, mask='@CA')
        traj_immut.superpose(ref=ref, mask='@CA')
        
        aa_eq(traj_mut.xyz, traj_immut.xyz)

        # test saving
        with tempfolder():
            traj_mut.save('t0.nc')
            traj_immut.save('t1.nc')

            aa_eq(pt.load('t0.nc', traj_mut.top).xyz,
                  pt.load('t1.nc', traj_immut.top).xyz)

        # turn off superpose
        traj_immut._is_superposed = False
        aa_eq(traj_immut.xyz, traj_immut2.xyz)

if __name__ == "__main__":
    unittest.main()
