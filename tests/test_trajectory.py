from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestAutoImageAndRotateDihedral(unittest.TestCase):

    def test_autoimage_rotatedihedral(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        farray = traj[:]

        t0trajectory = pt.Trajectory(traj)
        aa_eq(farray.unitcells, t0trajectory.unitcells)

        # autoimage
        farray.autoimage()
        t0trajectory.autoimage()
        aa_eq(farray.xyz, t0trajectory.xyz)

        # rotate_dihedral
        pt.rotate_dihedral(t0trajectory, '3:phi:120')
        pt.rotate_dihedral(farray, '3:phi:120')
        aa_eq(farray.xyz, t0trajectory.xyz)

        aa_eq(
            pt.calc_phi(t0trajectory, '3').values,
            [120 for _ in range(t0trajectory.n_frames)])


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
        saved_traj = pt.iterload("./data/tz2.autoimage.nc",
                                 "./data/tz2.ortho.parm7")

        # make sure to reproduce cpptraj's output too
        aa_eq(saved_traj.xyz, fa.xyz)

    def test_FromIterable(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        aa_eq(pt.Trajectory.from_iterable(traj).xyz, traj.xyz)


class TestAppend(unittest.TestCase):

    def test_append_trajectory(self):
        # test append
        traj = pt.Trajectory()
        t = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        traj.top = t.top

        # append single Frame
        traj.append(t[0])
        assert traj.n_frames == 1

        # append xyz
        traj.append(t.xyz[:])
        assert traj.n_frames == t.n_frames + 1

        # append TrajectoryIterator
        traj.append(t)
        assert traj.n_frames == t.n_frames * 2 + 1

        # append frame_iter
        traj.append(t.iterframe())
        assert traj.n_frames == t.n_frames * 3 + 1

        # append pt.iterframe_master
        traj.append(pt.iterframe_master(t))
        assert traj.n_frames == t.n_frames * 4 + 1

        # append itself
        NFrames = traj.n_frames
        traj.append(traj)
        assert traj.n_frames == NFrames * 2

        # append itself frame_iter
        traj.append(traj.iterframe(stop=2))
        assert traj.n_frames == NFrames * 2 + 2

        # append pt.iterframe_master for itself
        NFrames = traj.n_frames
        traj.append(pt.iterframe_master(traj))
        assert traj.n_frames == NFrames * 2

        # append pt.iterframe_master for itself + other
        n0 = traj.n_frames
        n1 = t.n_frames
        traj.append(pt.iterframe_master([traj, t]))
        assert traj.n_frames == 2 * n0 + n1


class TestTrajectory(unittest.TestCase):

    def test_raise_construtor(self):
        self.assertRaises(ValueError, lambda: pt.Trajectory(pt.trajectory))
        # raise if filename is not string or list of string
        empty_traj = pt.Trajectory()
        empty_traj.top = pt.datafiles.load_tz2()[:].top
        self.assertRaises(ValueError, lambda: empty_traj.load(pt.Trajectory))

        # raise if empty Topology
        xyz = np.arange(90).astype('f8').reshape(3, 10, 3)
        self.assertRaises(ValueError, lambda: pt.Trajectory(xyz=xyz, top=pt.Topology()))

    def test_slice_basic(self):
        traj2 = pt.Trajectory()
        traj2.top = pt.load_topology("./data/Tc5b.top")
        traj2.load("./data/Tc5b.x")
        traj2.load("./data/Tc5b.x")
        traj2.load("./data/Tc5b.x")
        traj2.load("./data/Tc5b.x")
        fsub = traj2[2:10]
        fsub[0][0] = 100.

    def test_indexing(self):
        traj = pt.iterload('data/Tc5b.x', 'data/Tc5b.top')
        traj2 = pt.TrajectoryIterator()
        traj2.top = pt.load_topology("./data/Tc5b.top")
        traj2._load("./data/Tc5b.x")
        farray = traj2[[0, 9, 1]]
        assert farray.n_frames == 3
        assert traj2[0].atom(0) == farray[0].atom(0)
        assert traj2[9].atom(0) == farray[1].atom(0)
        assert traj2[1].atom(0) == farray[2].atom(0)

        arr = np.asarray(traj2[0]._buffer1d[:])
        frame0 = traj2[0]
        arr0 = np.asarray(frame0._buffer1d[:])

        mat0 = np.asmatrix(arr0).reshape(304, 3)
        mat0[:, 0] = np.asmatrix(list(range(304))).reshape(304, 1)
        assert frame0[0, 0] == 0.
        assert frame0[1, 0] == 1.
        assert frame0[2, 0] == 2.

        # raise if size = 0
        traj3 = pt.Trajectory()
        assert traj3.n_frames == 0, 'empty Trajectory, n_frames must be 0'
        self.assertRaises(IndexError, lambda: traj3[0])
        self.assertRaises(IndexError, lambda: traj3.__setitem__(0, traj[3]))

    def test_indexing_at(self):
        traj = pt.iterload('data/Tc5b.x', 'data/Tc5b.top')
        t0 = traj[:]

        aa_eq(t0.at(3).xyz, t0[3].xyz)
        aa_eq(t0.at(slice(None, None, None)).xyz, t0[slice(None, None, None)].xyz)

        at_3 = traj.at(3)
        aa_eq(at_3.xyz, traj[3].xyz)

        at_all = traj.at(slice(None, None, None))
        aa_eq(at_all.xyz, traj[slice(None, None, None)].xyz)

    def test_iter_basic(self):
        traj = pt.TrajectoryIterator()
        traj.top = pt.load_topology("./data/Tc5b.top")
        traj._load("./data/Tc5b.x")
        for frame in traj:
            pass

    def test_traj_topology(self):
        traj = pt.TrajectoryIterator()
        assert traj.top.is_empty() == True
        traj.top = pt.load_topology("./data/Tc5b.top")
        assert traj.top.is_empty() == False
        traj._load("./data/Tc5b.x")

        # use toplogy
        traj = pt.TrajectoryIterator()
        assert traj.topology.is_empty() == True
        traj.topology = pt.load_topology("./data/Tc5b.top")
        assert traj.topology.is_empty() == False
        traj._load("./data/Tc5b.x")

    def test_xyz(self):
        traj = pt.datafiles.load_tz2_ortho()
        t0 = traj[:]

        def set_xyz_not_c_contiguous():
            t0.xyz = np.asfortranarray(traj.xyz)

        def append_2d():
            traj1 = pt.load_sample_data('ala3')

        def set_xyz_not_same_n_atoms():
            traj1 = pt.load_sample_data('ala3')
            t0.xyz = traj1.xyz

        def append_2d():
            traj1 = pt.load_sample_data('ala3')
            t0.append_xyz(pt.tools.as_2darray(traj))

        # fortran order, autoconvert
        # make sure there is no TypeError
        set_xyz_not_c_contiguous()

        self.assertRaises(ValueError, lambda: set_xyz_not_same_n_atoms())
        self.assertRaises(ValueError, lambda: append_2d())

        # make sure to autoconvert from f4 to f8
        xyz_f4 = np.array(traj.xyz, dtype='f4')
        assert xyz_f4.itemsize == 4, 'must be f4'
        t0.xyz = xyz_f4
        aa_eq(t0.xyz, xyz_f4)
        assert t0.xyz.itemsize == 8, 'must be converted from f4 to f8'

    def test_from_iterables(self):
        '''test_from_iterables, tests are ind its doc. Test raise here
        '''
        traj = pt.datafiles.load_tz2_ortho()
        fi = pt.pipe(traj, ['autoimage'])
        # does not have Topology info
        self.assertRaises(ValueError, lambda: pt.Trajectory.from_iterable(fi))

    def test_add_merge_two_trajs(self):
        '''test_add_merge_two_trajs'''
        traj1 = pt.datafiles.load_ala3()[:]
        traj2 = pt.datafiles.load_rna()[:]
        trajiter = pt.datafiles.load_rna()(stop=1)

        # raise if do not have the same n_frames
        self.assertRaises(ValueError, lambda: traj1 + traj2)

        traj1 = traj1[:1]
        traj2 = traj2[:1]

        traj3 = traj1 + traj2

        aa_eq(traj1.xyz, traj3.xyz[:, :traj1.n_atoms])
        aa_eq(traj2.xyz, traj3.xyz[:, traj1.n_atoms:])
        aa_eq(pt.tools.merge_trajs(traj1, traj2).xyz, traj3.xyz)
        aa_eq(pt.tools.merge_trajs(traj1, traj2).xyz, (traj1 + trajiter).xyz)

    def test_allocate_frames(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        traj2 = pt.Trajectory()
        traj2._allocate(traj.n_frames, traj.n_atoms)
        assert (traj2.shape == traj.shape)

        traj2.top = traj.top.copy()
        traj2.xyz = traj.xyz[:]
        aa_eq(traj2.xyz, traj.xyz)


class TestSaveToDisk(unittest.TestCase):

    def test_basic_saving(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        fa = traj[:]
        fname = "./output/test_savemethod.x"
        fname2 = "./output/test_savemethod_2.x"
        fa.save(fname, overwrite=True)
        traj.save(fname2, overwrite=True)

        # load
        fanew = pt.iterload(fname, fa.top)
        fanew2 = pt.iterload(fname2, fa.top)
        assert fanew.n_frames == fa.n_frames == fanew2.n_frames

        for idx, f0 in enumerate(fa):
            f0new = fanew[idx]
            f0new2 = fanew2[idx]
            aa_eq(f0.xyz, f0new.xyz)
            aa_eq(f0.xyz, f0new2.xyz)

    def test_fancy_save(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        traj[1:8].save("./output/test_fancy_save_frame1_7.x", overwrite=True)

        fanew = pt.iterload("./output/test_fancy_save_frame1_7.x", traj.top)

        for idx, f0 in enumerate(traj[1:8]):
            f0new = fanew[idx]
            aa_eq(f0.xyz, f0new.xyz)


class TestSetitem(unittest.TestCase):

    def test_setitem(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]

        # single value
        fa[0, 0, 0] = 100.
        assert fa[0, 0, 0] == 100.

        # single atom
        fa[0, 0] = [100., 101, 102]
        aa_eq(fa[0, 0], [100, 101, 102])

        # a set of atoms
        indices = [1, 10, 11]
        mask = '@2,11,12'
        xyz_sub = fa.xyz[:, indices] + 1.
        fa[mask] = xyz_sub
        aa_eq(fa[mask].xyz, xyz_sub)

        # all atoms
        xyz = traj.xyz + 2.
        fa["*"] = xyz
        aa_eq(fa.xyz, xyz)

        # all atoms for a set of frames
        xyz = traj.xyz[:3] + 3.
        fa[:3]["*"] = xyz
        aa_eq(fa.xyz[:3], xyz)

        # automatically cast
        fa0 = fa.copy()
        xyz = fa.xyz + 1.
        fa0[0] = xyz[0]  # fa[0] return a Frame
        aa_eq(fa0[0].xyz, xyz[0])
        # try to assign a Frame
        fa0[0] = fa[0]
        aa_eq(fa0[0].xyz, fa[0].xyz)

        def shape_mismatch():
            fa[0] = xyz

        self.assertRaises(ValueError, lambda: shape_mismatch())

        def shape_mismatch2():
            fa[0] = pt.Frame()

        self.assertRaises(ValueError, lambda: shape_mismatch2())

        # assign to None
        def None_value():
            fa[0] = None

        self.assertRaises(ValueError, lambda: None_value())


if __name__ == "__main__":
    unittest.main()
