#!/usr/bin/env python
import unittest
import pytraj as pt
from pytraj import Topology, Trajectory, TrajectoryIterator
from pytraj.testing import aa_eq



class TestIO(unittest.TestCase):
    def setUp(self):
        self.traj_tz2_ortho = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

    def test_save_traj_from_file(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:5]
        pt.write_traj(filename="./output/test_0.binpos",
                        traj=traj,
                        top="./data/Tc5b.top",
                        overwrite=True)

        savedtraj = pt.iterload("./output/test_0.binpos", traj.top)
        assert savedtraj.n_frames == traj.n_frames

    def test_blindload(self):
        top = pt.load_topology("./data/Tc5b.top")
        assert isinstance(top, Topology) == True

        traj = pt.iterload(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")

        is_traj = (
            isinstance(traj, TrajectoryIterator) or
            isinstance(traj, Trajectory))
        assert is_traj == True

    def test_ParmFile(self):
        top = pt.read_parm("./data/Tc5b.top")
        pt.write_parm("./output/test_io.top", top)
        newtop = pt.read_parm("./output/test_io.top")
        assert top.n_atoms == newtop.n_atoms

    def test_load_and_save_0(self):
        # need to load to Trajectory to save
        traj = pt.iterload(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")[:]

        indices = list(range(2, 3, 5)) + [3, 7, 9, 8]
        pt.write_traj(filename="./output/test_io_saved_.x",
                        traj=traj[:],
                        top="./data/Tc5b.top",
                        frame_indices=indices,
                        overwrite=True)

        # check frames
        traj2 = pt.iterload(
            filename="./output/test_io_saved_.x",
            top="./data/Tc5b.top")

        # about 50% failures
        assert traj2.n_frames == len(indices)

    def test_load_and_save_1(self):
        traj = pt.iterload(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")

        indices = list(range(2, 4)) + [3, 7, 9, 8]
        pt.write_traj(filename="./output/test_io_saved.pdb",
                        traj=traj,
                        top="./data/Tc5b.top",
                        frame_indices=indices,
                        overwrite=True)

        # check frames
        traj = pt.iterload(
            filename="./output/test_io_saved.pdb",
            top="./data/Tc5b.top")
        assert traj.n_frames == len(indices)
        assert traj.top.n_atoms == 304

    def test_get_coordinates_trajecotoryiterator(self):
        '''immutable pytraj.TrajectoryIterator
        '''
        traj = self.traj_tz2_ortho.copy()

        # all coordinates
        xyz = pt.get_coordinates(traj)
        aa_eq(traj.xyz, xyz)

        # given frames
        xyz = pt.get_coordinates(traj, frame_indices=[0, 5])
        aa_eq(traj[[0, 5]].xyz, xyz)

        # given frames, autoimage=True
        xyz = pt.get_coordinates(traj, frame_indices=[0, 5], autoimage=True)
        aa_eq(traj[[0, 5]].autoimage().xyz, xyz)

        # given frames, autoimage=True, rmsfit=ref
        ref = traj[-3]
        pt.autoimage(ref, top=traj.top)
        xyz = pt.get_coordinates(traj, frame_indices=[0, 5], autoimage=True, rmsfit=ref)
        aa_eq(traj[[0, 5]].autoimage().superpose(ref).xyz, xyz)

        xyz = pt.get_coordinates(traj, frame_indices=range(3), autoimage=True, rmsfit=ref)

    def test_get_coordinates_trajecotory(self):
        '''mutable pytraj.Trajectory
        '''
        traj = self.traj_tz2_ortho.copy()
        # make a different copy since ``traj`` is mutable
        traj2 = traj.copy()

        # all coordinates
        xyz = pt.get_coordinates(traj)
        aa_eq(traj.xyz, xyz)

        # given frames
        xyz = pt.get_coordinates(traj, frame_indices=[0, 5])
        aa_eq(traj[[0, 5]].xyz, xyz)

        # given frames, autoimage=True
        xyz = pt.get_coordinates(traj, frame_indices=[0, 5], autoimage=True)
        aa_eq(traj2[[0, 5]].autoimage().xyz, xyz)

        # given frames, autoimage=True, rmsfit=ref
        ref = traj[-3]
        pt.autoimage(ref, top=traj.top)
        xyz = pt.get_coordinates(traj, frame_indices=[0, 5], autoimage=True, rmsfit=ref)
        aa_eq(traj2[[0, 5]].autoimage().superpose(ref).xyz, xyz)


if __name__ == "__main__":
    unittest.main()
