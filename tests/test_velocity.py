#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import nose.tools as nt


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

    def test_velocityautocorr(self):
        # usevelocity = False
        traj = pt.iterload("./data/issue807/trunc.nc",
                           "data/issue807/system.prmtop")
        data0 = pt.all_actions.velocityautocorr(traj, tstep=2, norm=True, direct=True)

        cm = """
        parm data/issue807/system.prmtop
        trajin data/issue807/trunc.nc
        velocityautocorr mydata * tstep 2 norm direct
        """
        
        state = pt.load_cpptraj_state(cm)
        state.run()
        aa_eq(data0, state.data[-1].values)

        # usevelocity = True
        data = pt.all_actions.velocityautocorr(traj, 
                tstep=1, norm=False, direct=True, usevelocity=True)

        cm = """
        parm data/issue807/system.prmtop
        trajin data/issue807/trunc.nc
        velocityautocorr mydata * tstep 1 direct usevelocity
        """
        
        state = pt.load_cpptraj_state(cm)
        state.run()
        aa_eq(data, state.data[-1].values)

        # try on memory Trajectory, usevelocity = False
        traj_on_mem = traj[:]
        data2 = pt.all_actions.velocityautocorr(traj_on_mem, tstep=2, norm=True, direct=True)
        aa_eq(data0, data2)

        # try on memory Trajectory, usevelocity = True
        # need to raise if no `velocity_arr` is given
        traj_on_mem = traj[:]

        def func():
            pt.all_actions.velocityautocorr(traj_on_mem, tstep=2, norm=True, direct=True,
                                     usevelocity=True)
        nt.assert_raises(ValueError, func)

        # try on memory Trajectory, usevelocity = True
        traj_on_disk = traj
        velocity_arr = pt.get_velocity(traj_on_disk)
        data3 = pt.all_actions.velocityautocorr(traj_on_mem, tstep=2, norm=True, direct=True,
                                     usevelocity=True,
                                     velocity_arr=velocity_arr)

        data4 = pt.all_actions.velocityautocorr(traj_on_disk, tstep=2, norm=True, direct=True,
                                     usevelocity=True)

        # raise if velocity_arr has wrong shape
        velocity_arr_2 = velocity_arr.flatten()
        def func2():
            pt.all_actions.velocityautocorr(traj_on_mem, usevelocity=True, velocity_arr=velocity_arr_2)
        nt.assert_raises(ValueError, func2)

if __name__ == "__main__":
    unittest.main()
