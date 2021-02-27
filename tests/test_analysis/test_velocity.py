#!/usr/bin/env python
from __future__ import print_function
import os
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.testing import cpptraj_test_dir
import pytest


class TestVelocity(unittest.TestCase):
    def test_set_velocity(self):
        traj = pt.load(
            os.path.join(cpptraj_test_dir, 'tz2.rst7'), fn('tz2.parm7'))
        saved_rst7 = pt.iterload(
            os.path.join(cpptraj_test_dir, 'Test_SetVelocity',
                         'tz2.vel.rst7.save'), fn('tz2.parm7'))
        # Set default RNG back to Marsaglia
        pt.set_default_rng(0)
        pt.set_velocity(traj, temperature=298, ig=10)
        aa_eq(traj[0].velocity, saved_rst7[0].velocity)

    def test_get_velocity(self):
        traj = pt.iterload(
            fn('issue807/trunc.nc'), fn("issue807/system.prmtop"))

        traj[0]

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
        traj2 = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
        with pytest.raises(ValueError):
            pt.get_velocity(traj2)

    def test_velocityautocorr(self):
        # usecoords = False
        traj = pt.iterload(
            fn('issue807/trunc.nc'),
            fn("issue807/system.prmtop"),
            frame_slice=(0, 3))
        data0 = pt.all_actions.velocityautocorr(
            traj, tstep=2, norm=True, direct=True)

        cm = """
        parm {}
        trajin {} 1 3
        velocityautocorr mydata * tstep 2 norm direct
        """.format(fn('issue807/system.prmtop'), fn('issue807/trunc.nc'))

        state = pt.load_cpptraj_state(cm)
        state.run()
        aa_eq(data0[0], state.data[-2].values)

        # usecoords = True
        data = pt.all_actions.velocityautocorr(
            traj, tstep=1, norm=False, direct=True, usecoords=True)

        cm = """
        parm {}
        trajin {} 1 3
        velocityautocorr mydata * tstep 1 direct usecoords
        """.format(fn('issue807/system.prmtop'), fn('issue807/trunc.nc'))

        state = pt.load_cpptraj_state(cm)
        state.run()
        aa_eq(data[0], state.data[-2].values)

        # try on memory Trajectory, usecoords = False
        traj_on_mem = traj[:]
        data2 = pt.all_actions.velocityautocorr(
            traj_on_mem, tstep=2, norm=True, direct=True)
        aa_eq(data0[0], data2[0])

        # try on memory Trajectory, usecoords = True
        # need to raise if no `velocity_arr` is given
        traj_on_mem = traj[:]

        # try on memory Trajectory, usecoords = True
        traj_on_disk = traj
        velocity_arr = pt.get_velocity(traj_on_disk)
        data3 = pt.all_actions.velocityautocorr(
            traj_on_mem,
            tstep=2,
            norm=True,
            direct=True,
            usecoords=True,
            velocity_arr=velocity_arr)

        data4 = pt.all_actions.velocityautocorr(
            traj_on_disk, tstep=2, norm=True, direct=True, usecoords=True)

        # raise if velocity_arr has wrong shape
        velocity_arr_2 = velocity_arr.flatten()
        with pytest.raises(ValueError):
            pt.all_actions.velocityautocorr(
                traj_on_mem, usecoords=True, velocity_arr=velocity_arr_2)
