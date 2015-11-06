#!/usr/bin/env python
import unittest
import numpy as np
import pytraj as pt
from pytraj import Topology, Trajectory, TrajectoryIterator
from pytraj.testing import aa_eq


try:
    import scipy
    has_scipy = True
except ImportError:
    has_scipy = False

class TestIO(unittest.TestCase):
    def setUp(self):
        self.traj_tz2_ortho = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

    def test_load_comprehensive(self):
        traj = self.traj_tz2_ortho
        fn, tn = ("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

        # frame_slice
        t0 = pt.io.load_traj(fn, tn, frame_slice=(0, 3))
        aa_eq(self.traj_tz2_ortho[:3].xyz, t0.xyz)

        # mask
        t1 = pt.load(fn, tn, mask='@CA')
        aa_eq(t1.xyz, traj['@CA'].xyz)

        # frame_indices, list
        t1 = pt.load(fn, tn, frame_indices=[0, 3])
        aa_eq(t1.xyz, traj[[0, 3]].xyz)

        # frame_indices, tuple
        t1 = pt.load(fn, tn, frame_indices=(0, 3))
        aa_eq(t1.xyz, traj[[0, 3]].xyz)

        # mask and frame_indices
        t2 = pt.load(fn, tn, mask='@CA', frame_indices=[3, 8])
        aa_eq(t2.xyz, traj[[3, 8], '@CA'].xyz)

        # load http
        t2 = pt.load('https://raw.githubusercontent.com/ParmEd/ParmEd/master/test/files/2koc.pdb')
        assert t2.n_atoms == 451, '2koc, 451 atoms'
        assert isinstance(t2, pt.Trajectory), 'must be Trajectory when loading 2koc.pdb'


    def test_save_traj_from_file(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:5]
        pt.write_traj(filename="./output/test_0.binpos",
                        traj=traj,
                        top="./data/Tc5b.top",
                        overwrite=True)

        savedtraj = pt.iterload("./output/test_0.binpos", traj.top)
        assert savedtraj.n_frames == traj.n_frames

        # write_xyz
        pt.write_traj("./output/test_0.nc",
                      traj.xyz, 
                      top="./data/Tc5b.top",
                      overwrite=True)
        aa_eq(pt.iterload('output/test_0.nc', traj.top).xyz, traj.xyz)

        # write single Frame
        pt.write_traj("./output/test_0.nc",
                      traj[0],
                      top=traj.top,
                      overwrite=True)
        aa_eq(pt.iterload('output/test_0.nc', traj.top).xyz, traj[0].xyz)

        # raise if traj is None
        self.assertRaises(ValueError, lambda: pt.write_traj("./output/test_0.nc",
                      None,
                      overwrite=True))

        # raise if _top is None
        fi = pt.create_pipeline(traj, ['autoimage', ])
        self.assertRaises(ValueError, lambda: pt.write_traj("./output/test_0.nc",
                      traj=fi,
                      overwrite=True))

        # raise if Frame with frame_indices
        self.assertRaises(ValueError, lambda: pt.write_traj("./output/test_0.nc",
                      traj[0],
                      top="./data/Tc5b.top",
                      frame_indices=[3, 2],
                      overwrite=True))

        # raise if Frame with no Topology
        self.assertRaises(ValueError, lambda: pt.write_traj("./output/test_0.nc",
                      traj[0],
                      overwrite=True))

        # test if xyz is not c-contiguous
        # pytraj will autoconvert to c-contiguous
        xyz = np.asfortranarray(traj.xyz)
        # make sure no ValueError or TypeError is raised
        pt.write_traj('output/xyz.nc', xyz, top=traj.top, overwrite=True)

    def test_blind_load(self):
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
        # do not support frame_indices for TrajectoryIterator
        self.assertRaises(ValueError, lambda: pt.iterload(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top", frame_indices=[0, 5]))

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

        # with mask
        xyz = pt.get_coordinates(traj, mask='@CA')
        aa_eq(xyz, traj['@CA'].xyz)

        # slow
        fi = pt.create_pipeline(traj, ['autoimage'])
        aa_eq(pt.get_coordinates(fi), traj[:].autoimage().xyz)

        # raise
        self.assertRaises(ValueError, lambda: pt.get_coordinates(traj(), frame_indices=[0,
            2]))

    def test_get_coordinates_trajectory(self):
        '''mutable pytraj.Trajectory
        '''
        traj = pt.Trajectory(xyz=self.traj_tz2_ortho.xyz, top=self.traj_tz2_ortho.top)
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

    def test_load_and_save_topology(self):
        top = self.traj_tz2_ortho.top
        top.save('output/test.prmtop')
        top2 = pt.load_topology('output/test.prmtop')
        assert top.n_atoms == top2.n_atoms, 'must have the same n_atoms'

        # shortcut
        pt.load_topology('data/tz2.parm7').save('output/tz2_0.parm7')
        top3 = pt.load_topology('output/tz2_0.parm7')
        assert top3.n_atoms == 223, 'must have 223 atoms'
        assert top3.n_residues == 13, 'must have 13 residues'

        # raise
        self.assertRaises(ValueError, lambda: pt.load_topology(100))

    @unittest.skipIf(not has_scipy, 'skip since does not have scipy')
    def test_load_netcdf(self):
        fname = 'data/tz2.ortho.nc'
        tname = 'data/tz2.ortho.parm7'
        saved_traj = pt.iterload(fname, tname)

        # load all
        traj = pt.io._load_netcdf(fname, tname)
        aa_eq(saved_traj.xyz, traj.xyz)
        aa_eq(saved_traj.unitcells, traj.unitcells)

        # load some
        traj = pt.io._load_netcdf(fname, tname, frame_indices=range(5))
        aa_eq(saved_traj[range(5)].xyz, traj.xyz)
        aa_eq(saved_traj.unitcells[range(5)], traj.unitcells)

    def test_single_frame(self):
        traj = pt.load_sample_data('tz2')
        frame = pt.io.load_frame(traj.filename, traj.top.filename, 3)
        aa_eq(traj[3].xyz, frame.xyz)

        self.assertRaises(RuntimeError, lambda: pt.io.load_frame('afddsfdsfa',
            traj.top.filename, 3))

        self.assertRaises(RuntimeError, lambda: pt.io.load_frame(filename='afddsfdsfa',
            top=traj.top.filename, index=3))

    def test_download_pdb(self):
        pt.io.download_PDB('1l2y', 'output/', overwrite=True)
        t2 = pt.load('output/1l2y.pdb')
        assert t2.n_atoms == 304, 'must have 304 atoms'
        self.assertRaises(ValueError, lambda: pt.io.download_PDB('1l2y', 'output/',
            overwrite=False))

    def test_load_https(self):
        top = pt.io.load_topology('https://raw.githubusercontent.com/ParmEd/ParmEd/master/test/files/2koc.pdb')
        assert top.n_atoms == 451, '2koc'


if __name__ == "__main__":
    unittest.main()