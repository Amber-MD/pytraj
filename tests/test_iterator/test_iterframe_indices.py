from __future__ import print_function
import pytraj as pt
from utils import fn
import numpy as np
from pytraj.testing import aa_eq
import pytest


class TestIteraframeIndices:
    def test_iterframe_indices(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        t0 = traj[:]
        indices = range(3)

        d1 = pt.radgyr(traj[indices])
        d2 = pt.radgyr(traj, frame_indices=indices)

        aa_eq(d2, d1)

        n_frames = traj.n_frames
        aa_eq(
             np.array([frame.xyz.copy() for frame in
                       traj._iterframe_indices([-2, -1])]),
             traj[[n_frames-2, n_frames-1]].xyz) # yapf: disable
        aa_eq(
             traj[[-2, -1]].xyz,
             traj[[n_frames-2, n_frames-1]].xyz) # yapf: disable
