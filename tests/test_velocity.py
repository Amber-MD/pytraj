#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestVelocity(unittest.TestCase):

    def test_velocity(self):
        traj = pt.iterload("./data/issue807/trunc.nc",
                           "data/issue807/system.prmtop")

        f = traj[0]

        # no mask, no frame_indices
        vels = pt.get_velocity(traj)
        assert vels.shape == (traj.n_frames, traj.n_atoms, 3), 'vels.shape'

        # string mask
        vels = pt.get_velocity(traj, '@O', frame_indices=[0, 2])
        fi = traj(frame_indices=[0, 2], mask='@O')
        assert vels.shape == (fi.n_frames, fi.top.n_atoms, 3), 'vels.shape'

        # atom indices
        atm_indices = pt.select_atoms('@O', traj.top)
        vels_ = pt.get_velocity(traj, atm_indices, frame_indices=[0, 2])
        fi = traj(frame_indices=[0, 2], mask='@O')
        assert vels_.shape == (fi.n_frames, fi.top.n_atoms, 3), 'vels.shape'
        aa_eq(vels, vels_)

        # raise if not having velocity
        traj2 = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
        self.assertRaises(ValueError, lambda: pt.get_velocity(traj2))


if __name__ == "__main__":
    unittest.main()
