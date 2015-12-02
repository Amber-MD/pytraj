#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestREMDTemperature(unittest.TestCase):

    def setUp(self):
        self.trajin_text = '''
        parm  data/Test_RemdTraj/ala2.99sb.mbondi2.parm7
        trajin data/Test_RemdTraj/rem.nc.000 remdtraj remdtrajtemp 300.
        distance @10 @20
        '''

    def test_load_cpptraj_state_from_text(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        text = '''
        parm  data/Test_RemdTraj/ala2.99sb.mbondi2.parm7
        trajin data/Test_RemdTraj/rem.nc.000 remdtraj remdtrajtemp 300.
        distance @10 @20
        '''

        input_file = 'data/Test_RemdTraj/traj.in'

        state_from_file = pt.load_cpptraj_state(input_file)
        state_from_file.run()
        # remove DatasetTopology
        data = state_from_file.data
        data.remove_set(data[0])

        data_0 = state_from_file.data.values

        state_from_text = pt.datafiles.load_cpptraj_state(text)
        state_from_text.run()
        data = state_from_text.data
        data.remove_set(data[0])
        data_1 = state_from_text.data.values
        aa_eq(data_0, data_1)

    def test_iterload_and_load_remd(self):
        # iterload_remd
        traj = pt.iterload_remd("data/Test_RemdTraj/rem.nc.000",
                                "data/Test_RemdTraj/ala2.99sb.mbondi2.parm7",
                                T=300.0)
        for frame in traj:
            assert frame.temperature == 300.0, 'frame temperature must be 300.0 K'
        dist = pt.distance(traj, '@10 @20')

        state = pt.load_cpptraj_state(self.trajin_text)
        state.run()
        aa_eq(dist, state.data[-1])

        # load_remd
        traj2 = pt.load_remd("data/Test_RemdTraj/rem.nc.000",
                             "data/Test_RemdTraj/ala2.99sb.mbondi2.parm7",
                             T=300.0)
        aa_eq(traj.xyz, traj2.xyz)

        # with Topology
        traj2 = pt.iterload_remd("data/Test_RemdTraj/rem.nc.000",
                                 top=traj.top,
                                 T=300.0)
        aa_eq(traj.xyz, traj2.xyz)


if __name__ == "__main__":
    unittest.main()
