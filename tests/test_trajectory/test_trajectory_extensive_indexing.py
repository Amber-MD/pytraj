import unittest
import pytraj as pt
from utils import fn

from pytraj.testing import aa_eq
from pytraj import *

refilename = fn('Tc5b.nat.crd')
trajin = """
"""

ts = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

# create Trajectory to store Frame
FRAMENUM = 999
FARRAY = ts[:FRAMENUM]


class TestTrajectory:
    def test_len_in_memory(self):
        N = 10
        farray = FARRAY[:N].copy()
        assert farray.n_frames == N
        old_xyz_5_10 = farray[5].xyz[:10]
        assert farray[:3].n_frames == 3
        assert farray[1:3].n_frames == 2
        assert farray[3:1].n_frames == 0
        assert farray[3:1:-1].n_frames == 2
        assert farray[-1:-3].n_frames == 0
        assert farray[-1:-3:-1].n_frames == 2
        aa_eq(farray[-1].xyz, farray[N - 1].xyz)

        # need to create a temp farray
        subfarray = farray[5:1:-1]
        aa_eq(subfarray[0].xyz, farray[5].xyz)
        aa_eq(old_xyz_5_10, farray[5].xyz[:10])

        f_last = farray[-3:-1][-1]

    def test_len_on_disk(self):
        N = ts.n_frames
        old_xyz_5_10 = ts[5].xyz[:10].copy()
        assert ts[:3].n_frames == 3
        assert ts[1:3].n_frames == 2
        assert ts[3:1].n_frames == 0
        assert ts[3:1:-1].n_frames == 2
        assert ts[-1:-3].n_frames == 0
        assert ts[-1:-3:-1].n_frames == 2
        # need to store xyz
        xyz = ts[-1].xyz.copy()
        aa_eq(xyz, ts[N - 1].xyz)

        # need to create a temp farray
        subfarray = ts[5:1:-1]
        aa_eq(subfarray[0].xyz, ts[5].xyz)
        aa_eq(old_xyz_5_10, ts[5].xyz[:10])

    def test_mask_indexing(self):
        # In-memory Trajectory
        traj = ts[:]
        assert traj["@CA"].shape == (traj.n_frames, traj.top("@CA").n_atoms, 3)
        assert traj[2:4]["@CA"].shape == (2, traj.top("@CA").n_atoms, 3)

        # On-disk Trajectory
        traj = ts
        assert traj["@CA"].shape == (traj.n_frames, traj.top("@CA").n_atoms, 3)
        assert traj[2:4]["@CA"].shape == (2, traj.top("@CA").n_atoms, 3)

    def test_negative_indexing(self, tmpdir):
        mt = pt.load(fn('tz2.nc'), fn('tz2.parm7'))  # mt = in-memory
        dt = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))  # mt = on-disk
        aa_eq(mt[-1].xyz, dt[-1].xyz)

        # one-frame trajectory
        with tmpdir.as_cwd():
            dt[:1].save("one.nc")
            dt2 = pt.iterload("one.nc", dt.top)
            mt2 = pt.load("one.nc", dt.top)
            assert dt2.n_frames == mt2.n_frames == 1
            aa_eq(mt2[0].xyz, mt2[-1].xyz)
            aa_eq(dt2[0].xyz.copy(), dt2[-1].xyz.copy()) # FIXME: wrong result if not using copy
            # weird memory stuff
            aa_eq(mt2[-1].xyz, dt2[-1].xyz)
