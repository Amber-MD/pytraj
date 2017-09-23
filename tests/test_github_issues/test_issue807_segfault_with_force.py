#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn


class TestSegmetationFault(unittest.TestCase):
    def test_issue807(self):
        # files are provided by Chris Lee
        traj = pt.iterload(
            fn('issue807/trunc.nc'), fn("issue807/system.prmtop"))

        traj[0]
        for frame in traj:
            pass
        traj[::2]
        pt.radgyr(traj[[0, 3]])
        pt.radgyr(traj, frame_indices=[0, 3])
        pt.radgyr(traj())
        traj[:3, '@O'].xyz
        pt.get_coordinates((traj.iterframe(mask='@O')))
        pt.radgyr(traj(mask='@O'))
        for c in pt.iterchunk(traj, 4):
            assert c[0].n_atoms == traj.top.n_atoms
