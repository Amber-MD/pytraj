#!/usr/bin/env python
import os
import sys
import unittest
import numpy as np
import pytest
import pytraj as pt
from pytraj import Topology, Trajectory, TrajectoryIterator, TrajectoryWriter
from pytraj.testing import aa_eq
from pytraj.testing import cpptraj_test_dir
from pytraj.testing import get_remd_fn
from pytraj.testing import tempfolder
from pytraj.io import _get_amberhome

# local
from utils import fn

try:
    has_scipy = True
except ImportError:
    has_scipy = False

amberhome = os.getenv('AMBERHOME', '')
tleap = amberhome + '/bin/tleap'
traj_tz2_ortho = pt.iterload(fn("tz2.ortho.nc"), fn("tz2.ortho.parm7"))
tc5b_trajin = fn('Tc5b.x')
tc5b_top = fn('Tc5b.top')


def test_iterload_comprehensive():
    trajin, tn = fn("tz2.ortho.nc"), fn("tz2.ortho.parm7")

    # frame_slice
    t0 = pt.iterload(trajin, tn, frame_slice=(0, -1, 0))
    aa_eq(traj_tz2_ortho.xyz, t0.xyz)

    t0 = pt.iterload(trajin, tn, frame_slice=(0, -1, 2))
    aa_eq(traj_tz2_ortho.xyz[::2], t0.xyz)

    # stride
    t0 = pt.iterload(trajin, tn, stride=2)
    aa_eq(traj_tz2_ortho.xyz[::2], t0.xyz)

    # stride, ignore frame_slice
    t0 = pt.iterload(trajin, tn, stride=2, frame_slice=(0, -1, 3))
    aa_eq(traj_tz2_ortho.xyz[::2], t0.xyz)

    # stride, load two files
    t0 = pt.iterload([trajin, trajin], tn, stride=2)
    xyz_2 = np.vstack((traj_tz2_ortho.xyz[::2], traj_tz2_ortho.xyz[::2]))
    aa_eq(xyz_2, t0.xyz)

    # stride, load two files, ignore frame_slice
    t0 = pt.iterload(
        [trajin, trajin], tn, stride=2, frame_slice=[(0, -1, 5), (0, -1, 2)])
    xyz_2 = np.vstack((traj_tz2_ortho.xyz[::2], traj_tz2_ortho.xyz[::2]))
    aa_eq(xyz_2, t0.xyz)

    # stride, 4 trajs
    filenames, tn = get_remd_fn('remd_ala2')
    t0 = pt.iterload(filenames, tn, stride=3)
    # add frame_slice
    t1 = pt.iterload(
        filenames, tn, frame_slice=[
            (0, -1, 3),
        ] * 4)
    xyz_expected = np.vstack(
        [pt.iterload(fname, tn)[::3].xyz for fname in filenames])
    aa_eq(xyz_expected, t0.xyz)
    aa_eq(xyz_expected, t1.xyz)

    # stride, 4 trajs, ignore frame_slice
    filenames, tn = get_remd_fn('remd_ala2')
    t0 = pt.iterload(filenames, tn, stride=3, frame_slice=(0 - 1, 4))
    xyz_expected = np.vstack(
        [pt.iterload(fname, tn)[::3].xyz for fname in filenames])
    aa_eq(xyz_expected, t0.xyz)


def test_load_comprehensive():
    traj = traj_tz2_ortho
    trajin, tn = fn("tz2.ortho.nc"), fn("tz2.ortho.parm7")

    # load from filelist
    t0 = pt.load([trajin, trajin], tn)
    n_frames_half = int(t0.n_frames / 2)
    aa_eq(traj.xyz, t0[:n_frames_half].xyz)
    aa_eq(traj.xyz, t0[n_frames_half:].xyz)

    # frame_slice
    t0 = pt.io.load_traj(trajin, tn, frame_slice=(0, 3))
    aa_eq(traj_tz2_ortho[:3].xyz, t0.xyz)

    # mask
    t1 = pt.load(trajin, tn, mask='@CA')
    aa_eq(t1.xyz, traj['@CA'].xyz)

    # frame_indices, list
    t1 = pt.load(trajin, tn, frame_indices=[0, 3])
    aa_eq(t1.xyz, traj[[0, 3]].xyz)

    # frame_indices, tuple
    t1 = pt.load(trajin, tn, frame_indices=(0, 3))
    aa_eq(t1.xyz, traj[[0, 3]].xyz)

    # frame_indices with negative index
    t12 = pt.load(trajin, tn, frame_indices=(-2, -1))
    aa_eq(t12.xyz, traj[[traj.n_frames-2, traj.n_frames-1]].xyz)

    # mask and frame_indices
    t2 = pt.load(trajin, tn, mask='@CA', frame_indices=[3, 8])
    aa_eq(t2.xyz, traj[[3, 8], '@CA'].xyz)

    # stride
    t2 = pt.load(trajin, tn, stride=2)
    aa_eq(t2.xyz, traj[::2].xyz)

    # stride with mask
    t2 = pt.load(trajin, tn, stride=2, mask='@CA')
    aa_eq(t2.xyz, traj[::2, '@CA'].xyz)

    # stride, ignore frame_indices if stride is given
    t2 = pt.load(trajin, tn, stride=2, frame_indices=[2, 5, 8])
    aa_eq(t2.xyz, traj[::2].xyz)


def test_blind_load():
    top = pt.load_topology(tc5b_top)
    assert isinstance(top, Topology) == True

    traj = pt.iterload(filename=tc5b_trajin, top=tc5b_top)

    is_traj = (isinstance(traj, TrajectoryIterator)
               or isinstance(traj, Trajectory))
    assert is_traj


def test_ParmFile():
    top = pt.load_topology(tc5b_top)
    with tempfolder():
        pt.write_parm("test_io.top", top, overwrite=True)
        newtop = pt.load_topology("test_io.top")
        assert top.n_atoms == newtop.n_atoms
        # test raise if file exists
        with pytest.raises(RuntimeError):
            pt.write_parm("test_io.top", top, overwrite=False)


def test_load_and_save_0():
    # need to load to Trajectory to save
    traj = pt.iterload(filename=tc5b_trajin, top=tc5b_top)[:]

    indices = list(range(2, 3, 5)) + [3, 7, 9, 8]
    with tempfolder():
        pt.write_traj(
            filename="test_io_saved_.x",
            traj=traj[:],
            frame_indices=indices,
            overwrite=True)

        # check frames
        traj2 = pt.iterload(filename="test_io_saved_.x", top=tc5b_top)
        # about 50% failures
        assert traj2.n_frames == len(indices)


def test_load_and_save_1():
    # do not support frame_indices for TrajectoryIterator
    with pytest.raises(ValueError):
        pt.iterload(filename=tc5b_trajin, top=tc5b_top, frame_indices=[0, 5])

    traj = pt.iterload(filename=tc5b_trajin, top=tc5b_top)

    indices = list(range(2, 4)) + [3, 7, 9, 8]
    with tempfolder():
        pt.write_traj(
            filename="test_io_saved.pdb",
            traj=traj,
            frame_indices=indices,
            overwrite=True)

        # check frames
        traj = pt.iterload(filename="test_io_saved.pdb", top=tc5b_top)
        assert traj.n_frames == len(indices)
        assert traj.top.n_atoms == 304


def test_overwrite():
    trajin, tn = fn("tz2.nc"), fn("tz2.parm7")
    with tempfolder():
        traj = pt.iterload(trajin, tn)
        pt.write_traj('what.nc', traj)
        # ensure no IOError is raised.
        pt.write_traj('what.nc', traj, overwrite=True)


def test_get_coordinates_trajecotoryiterator():
    '''immutable pytraj.TrajectoryIterator
    '''
    traj = traj_tz2_ortho.copy()

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
    xyz = pt.get_coordinates(
        traj, frame_indices=[0, 5], autoimage=True, rmsfit=ref)
    aa_eq(traj[[0, 5]].autoimage().superpose(ref=ref).xyz, xyz)

    xyz = pt.get_coordinates(
        traj, frame_indices=range(3), autoimage=True, rmsfit=ref)

    # with mask
    xyz = pt.get_coordinates(traj, mask='@CA')
    aa_eq(xyz, traj['@CA'].xyz)

    # slow
    fi = pt.pipe(traj, ['autoimage'])
    aa_eq(pt.get_coordinates(fi), traj[:].autoimage().xyz)

    # raise
    with pytest.raises(ValueError):
        pt.get_coordinates(traj(), frame_indices=[0, 2])


def test_get_coordinates_trajectory():
    '''mutable pytraj.Trajectory
    '''
    traj = pt.Trajectory(xyz=traj_tz2_ortho.xyz, top=traj_tz2_ortho.top)
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
    xyz = pt.get_coordinates(
        traj, frame_indices=[0, 5], autoimage=True, rmsfit=ref)
    aa_eq(traj2[[0, 5]].autoimage().superpose(ref=ref).xyz, xyz)


def test_load_and_save_topology():
    top = traj_tz2_ortho.top
    with tempfolder():
        top.save('test.prmtop')
        top2 = pt.load_topology('test.prmtop')
        assert top.n_atoms == top2.n_atoms, 'must have the same n_atoms'

        # shortcut
        pt.load_topology(fn('tz2.parm7')).save('tz2_0.parm7')
        top3 = pt.load_topology('tz2_0.parm7')
        assert top3.n_atoms == 223, 'must have 223 atoms'
        assert top3.n_residues == 13, 'must have 13 residues'

        # raise
        with pytest.raises(ValueError):
            pt.load_topology(100)


def test_single_frame():
    traj = pt.load_sample_data('tz2')
    frame = pt.io.load_frame(traj.filename, traj.top.filename, 3)
    aa_eq(traj[3].xyz, frame.xyz)

    with pytest.raises(RuntimeError):
        pt.io.load_frame('afddsfdsfa', traj.top.filename, 3)

    with pytest.raises(RuntimeError):
        pt.io.load_frame(filename='afddsfdsfa', top=traj.top.filename, index=3)


# @unittest.skip('download_PDB')
def test_download_pdb():
    with tempfolder():
        pt.io.download_PDB('1l2y', './', overwrite=True)
        t2 = pt.load('1l2y.pdb')
        assert t2.n_atoms == 304, 'must have 304 atoms'
        with pytest.raises(ValueError):
            pt.io.download_PDB('1l2y', './', overwrite=False)


def test_short_save():
    with tempfolder():
        (pt.iterload(fn('tz2.nc'), fn('tz2.parm7')).save(
            'mini.nc', overwrite=True))
        assert pt.iterload(
            'mini.nc', fn('tz2.parm7')).n_frames == 101, 'must be 101 frames'


def test_formats():
    '''make sure to save correct format
    '''

    def assert_has_exptected_line_textfile(expected_line,
                                           line_index,
                                           fn='test'):
        with open(fn, 'r') as fh:
            lines = [fh.readline() for _ in range(3)]
        # print(lines)
        assert expected_line in lines[line_index], lines[line_index]

    def assert_has_exptected_line_binaryfile(expected_line, fn='test'):
        with open(fn, 'rb') as fh:
            buffer_ = fh.read()
        # print(buffer_)
        assert expected_line.encode() in buffer_

    traj = pt.datafiles.load_tz2()
    fn = 'test'

    with tempfolder():
        # pdb
        traj.save(fn, format='pdb', overwrite=True)
        expected_line = 'ATOM      1  N   SER'
        assert_has_exptected_line_textfile(expected_line, 0, fn)

        # mol2
        traj.save(fn, format='mol2', overwrite=True)
        expected_line = 'Cpptraj Generated mol2 file'
        assert_has_exptected_line_textfile(expected_line, 1, fn)

        # mdcrd
        traj.save(fn, format='mdcrd', overwrite=True)
        expected_line = 'Cpptraj Generated trajectory'
        assert_has_exptected_line_textfile(expected_line, 0, fn)

        # netcdf
        traj.save(fn, format='netcdf', overwrite=True)
        expected_line = 'CDF'
        assert_has_exptected_line_binaryfile(expected_line, fn)

        # dcd
        traj.save(fn, format='dcd', overwrite=True)
        expected_line = 'Cpptraj Generated dcd file'
        assert_has_exptected_line_binaryfile(expected_line, fn)

        # cif, not supported
        # traj.save(fn, format='cif', overwrite=True)
        # expected_line = 'Cpptraj Generated cif file'
        # assert_has_exptected_line_textfile(expected_line, 1, fn)

        # trr
        traj.save(fn, format='trr', overwrite=True)
        expected_line = 'GMX_trn_file'
        assert_has_exptected_line_binaryfile(expected_line, fn)


def test_load_url():
    url = 'https://raw.githubusercontent.com/Amber-MD/pytraj/master/tests/data/1L2Y.pdb'
    local_traj = pt.load(fn('1L2Y.pdb'))
    github_traj = pt.io.load_url(url)
    aa_eq(local_traj.xyz, github_traj.xyz)


def test_options():
    '''specify cpptraj options
    '''

    with tempfolder():
        # restart, keep extension
        traj_tz2_ortho[:1].save('test.rst7', options='keepext')
        assert os.path.exists('test.1.rst7')


def test_write_time():
    top_fname = os.path.join(cpptraj_test_dir, 'FtuFabI.NAD.TCL.parm7')
    traj_fname = os.path.join(cpptraj_test_dir, 'FtuFabI.NAD.TCL.nc')
    traj = pt.iterload(traj_fname, top_fname)
    time_arr = np.array([
        37031.0, 37032.0, 37033.0, 37034.0, 37035.0, 37036.0, 37037.0, 37038.0,
        37039.0, 37040.0
    ])
    aa_eq([f.time for f in traj], time_arr)
    traj2 = traj[:]  # in memory
    aa_eq(np.array([f.time for f in traj2]), time_arr)
    traj3 = traj2[:]  # slice in memory trajectory again
    aa_eq(np.array([f.time for f in traj3]), time_arr)

    with tempfolder():
        pt.write_traj('test.nc', traj, time=True)
        out_traj = pt.iterload('test.nc', traj.top)
        aa_eq(np.array([f.time for f in out_traj]), time_arr)

        # load method
        out_traj2 = pt.load('test.nc', traj.top)
        aa_eq(np.array([f.time for f in out_traj2]), time_arr)


def test_write_velocity_from_scratch():
    traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
    assert traj[0].velocity is None

    def add_velocity(traj):
        for frame in traj:
            frame.velocity = np.zeros(traj.n_atoms * 3) + 1.
            yield frame

    with tempfolder():
        out_fn = 'out.nc'
        with TrajectoryWriter(
                out_fn, top=traj.top, crdinfo={
                    'has_velocity': True
                }) as writer:
            for frame in add_velocity(traj):
                writer.write(frame)
        traj2 = pt.iterload(out_fn, top=traj.top)
        assert traj2.metadata['has_velocity']
        assert not traj2.metadata['has_force']
        # make sure no segmentation fault
        # issues/1486
        assert traj2[:].n_frames == traj2.n_frames


def test_write_force_from_scratch():
    traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
    assert traj[0].force is None

    def add_force(traj):
        for frame in traj:
            frame.force = np.zeros(traj.n_atoms * 3) + 1.
            yield frame

    with tempfolder():
        out_fn = 'out.nc'
        with TrajectoryWriter(
                out_fn, top=traj.top, crdinfo={
                    'has_force': True
                }) as writer:
            for frame in add_force(traj):
                writer.write(frame)
        traj2 = pt.iterload(out_fn, top=traj.top)
        assert traj2.metadata['has_force']
        assert not traj2.metadata['has_velocity']
        # issues/1486
        assert traj2[:].n_frames == traj2.n_frames


def test_write_both_force_and_velocity_from_scratch():
    traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
    assert traj[0].force is None
    assert traj[0].velocity is None

    def add_force_and_velocity(traj):
        for frame in traj:
            frame.force = np.zeros(traj.n_atoms * 3) + 1.
            frame.velocity = np.zeros(traj.n_atoms * 3) + 1.
            yield frame

    with tempfolder():
        out_fn = 'out.nc'
        with TrajectoryWriter(
                out_fn,
                top=traj.top,
                crdinfo={
                    'has_force': True,
                    'has_velocity': True
                }) as writer:
            for frame in add_force_and_velocity(traj):
                writer.write(frame)
        traj2 = pt.iterload(out_fn, top=traj.top)
        assert traj2.metadata['has_force']
        assert traj2.metadata['has_velocity']


def test_write_force_and_velocity():
    trajin = cpptraj_test_dir + '/Test_systemVF/systemVF.nc'
    tn = cpptraj_test_dir + '/Test_systemVF/systemVF.parm7'
    traj = pt.iterload(trajin, tn)
    assert traj.metadata['has_force']
    assert traj.metadata['has_velocity']

    with tempfolder():
        fn2 = 'test.nc'
        traj.save(fn2, overwrite=True, crdinfo=traj.metadata)

        traj2 = pt.iterload(fn2, traj.top)
        assert traj2.metadata['has_force']
        assert traj2.metadata['has_velocity']

        forces_traj = np.array([frame.force.copy() for frame in traj])
        forces_traj2 = np.array([frame.force.copy() for frame in traj2])
        aa_eq(forces_traj, forces_traj2)


def test_iterload_and_load_remd():
    # iterload_remd
    traj = pt.iterload_remd(
        fn("Test_RemdTraj/rem.nc.000"),
        fn("Test_RemdTraj/ala2.99sb.mbondi2.parm7"),
        T=300.0)
    for frame in traj:
        assert frame.temperature == 300.0, 'frame temperature must be 300.0 K'
    dist = pt.distance(traj, '@10 @20')

    trajin_text = '''
        parm  {}
        trajin {} remdtraj remdtrajtemp 300.
        distance @10 @20
    '''.format(
        fn('Test_RemdTraj/ala2.99sb.mbondi2.parm7'),
        fn('Test_RemdTraj/rem.nc.000 '))
    state = pt.load_cpptraj_state(trajin_text)
    state.run()
    aa_eq(dist, state.data[-1].values)

    # load_remd
    traj2 = pt.load_remd(
        fn("Test_RemdTraj/rem.nc.000"),
        fn("Test_RemdTraj/ala2.99sb.mbondi2.parm7"),
        T=300.0)
    aa_eq(traj.xyz, traj2.xyz)

    # with Topology
    traj2 = pt.iterload_remd(
        fn("Test_RemdTraj/rem.nc.000"), top=traj.top, T=300.0)
    aa_eq(traj.xyz, traj2.xyz)


def test_io_load_and_save_0():
    traj = pt.iterload(filename=tc5b_trajin, top=tc5b_top)[:10]
    indices = list(range(2, 3, 5)) + [3, 8, 9, 8]

    with tempfolder():
        pt.write_traj(
            filename="test_io_saved_.x",
            traj=traj,
            frame_indices=indices,
            overwrite=True)

        # check frames
        traj2 = pt.iterload(filename="test_io_saved_.x", top=tc5b_top)
        assert traj2.n_frames == len(indices)


def test_amberhome():
    if os.getenv('AMBERHOME') is None:
        fake_amberhome = '/home/hc/amber16'
        os.environ['AMBERHOME'] = fake_amberhome
        assert _get_amberhome() == fake_amberhome


@unittest.skipUnless(os.path.exists(tleap), 'must have tleap')
def test_leap():
    cm = """
    source leaprc.protein.ff14SB
    TC5b = sequence { NASN LEU TYR ILE GLN TRP LEU LYS
                            ASP GLY GLY PRO SER SER GLY ARG
                            PRO PRO PRO CSER }
    saveamberparm TC5b x.parm7 x.crd
    %s
    """

    for quit in ['quit', '']:
        verbose = False
        traj = pt.io.load_leap(cm % quit, verbose=verbose)
        assert traj.n_atoms == 304
