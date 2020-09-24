import pytest
import pytraj as pt
from utils import fn
import numpy as np
from pytraj import Trajectory
from pytraj.testing import aa_eq
from pytraj.utils import Timer
from pytraj.testing import cpptraj_test_dir
import time


class TestSlicingTrajectory:
    def test_array_like(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        xyz_source = traj.xyz[:].copy()
        traj_mem = traj[:]

        # slicing with list or array
        indices = [1, 2, 3]
        fa = traj[indices]
        fa2 = traj_mem[indices]
        fa3 = traj[range(1, 4)]
        fa4 = traj_mem[range(1, 4)]
        assert isinstance(fa, Trajectory)
        # from TrajectoryIterator
        aa_eq(fa[0].xyz, traj[1].xyz)
        aa_eq(fa[1].xyz, traj[2].xyz)
        # from Trajectory
        aa_eq(fa2[1].xyz, traj[2].xyz)
        aa_eq(fa2[0].xyz, traj[1].xyz)

        # a list with one element
        assert isinstance(traj[[1,]], Trajectory)
        aa_eq(traj[[1,]].xyz, traj[1:2].xyz)

        # from "range"
        aa_eq(fa3[1].xyz, traj[2].xyz)
        aa_eq(fa3[0].xyz, traj[1].xyz)
        aa_eq(fa4[1].xyz, traj[2].xyz)
        aa_eq(fa4[0].xyz, traj[1].xyz)

        # tuple
        aa_eq(traj[(1,)].xyz, xyz_source[1])

        # boolean
        bool_arr = np.array([False]*traj.n_frames)
        bool_arr[0] = True
        bool_arr[2] = True
        new_xyz = traj[bool_arr].xyz
        aa_eq(traj[bool_arr].xyz, traj_mem[bool_arr].xyz)
        aa_eq(new_xyz[0], traj[0].xyz)
        aa_eq(new_xyz[1], traj[2].xyz)

        with pytest.raises(IndexError):
            # Try to index with a list with len of 2 (< traj.n_frames)
            traj[[False, True]]


    def test_velocity(self):
        traj = pt.iterload(
            fn('issue807/trunc.nc'), fn("issue807/system.prmtop"))

        aa_eq(pt.get_velocity(traj), pt.get_velocity(traj[:]))

        aa_eq(
            pt.get_velocity(traj, frame_indices=[1, 3, 5]),
            pt.get_velocity(traj[[1, 3, 5]]))

    def test_atommask(self):
        # AtomMask
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        fa = traj[:]
        xyz = traj.xyz[:]
        atm = traj.top("@CA")
        indices = atm.indices

        aa_eq(fa[0, atm], fa[0][atm])
        aa_eq(traj[0, atm], fa[0][atm])
        aa_eq(traj[0, atm, 0], fa[0][atm, 0])
        aa_eq(traj[0, atm, 0], xyz[0][indices][0])


def test_slice_from_on_disk_trajectory():
    traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))[:]
    
    # list
    aa_eq(traj[3, 3], traj[3][3, :])

    # frame
    frame1 = traj[1]
    aa_eq(frame1[0], traj[1][:, :][0])
    assert traj[0, 0, 0] == -16.492
    assert traj[:, :, 0][0, 0] == traj[0, 0, 0]
    traj[0]
    farr0 = traj[:2]

    fa = traj[2:4]

    # we don't support traj[:, idx] or traj[:, idx, idy] since this give wrong answer
    #  got ~0.0 value
    aa_eq(traj[:, 0, 0], np.asarray(traj[0][0]))

    for i in range(traj[0]._buffer2d.shape[0]):
        aa_eq(traj[:, :, 0][i], traj[0]._buffer2d[i])

    # slicing with mask
    atm = traj.top("@CA")
    traj[atm]
    traj[:, atm]


def test_speed():

    top_fname = f"{cpptraj_test_dir}/DOPC.parm7"
    traj_fname = f"{cpptraj_test_dir}/DOPC.rst7"
    traj = pt.iterload([traj_fname]*100, top_fname)
    indices = list(range(50, 60))
    with Timer() as t0:
        xyz0 = traj[':WAT,OL', indices].xyz

    with Timer() as t1:
        xyz1 = traj[indices, ':WAT,OL'].xyz

    # https://github.com/Amber-MD/pytraj/issues/1494
    aa_eq(xyz0, xyz1)
    assert xyz0.shape[0] == 10
    assert abs(t0.value - t1.value) < 0.5


def test_segmentation_fault():
    # NOTE: no assert, just check for segfault
    traj = pt.load(fn('Tc5b.x'), fn('Tc5b.top'))
    pt.load(fn('Tc5b.x'), fn('Tc5b.top'))
    traj.top("@CA")
    f0 = traj[5]
    f0 = traj[0]
    f0.top = traj.top
    f0['@CA']
    traj[0, '@CA']

    f0 = traj[0, '@CA']
    f1 = traj['@CA'][0]
    assert pt.tools.rmsd(f0.xyz, f1.xyz) == 0.0
