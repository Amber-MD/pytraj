#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestREMDTemperature(unittest.TestCase):
    def test_load_cpptraj_state_from_text(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        text = '''
        parm  data/Test_RemdTraj/ala2.99sb.mbondi2.parm7
        trajin data/Test_RemdTraj/rem.nc.000 remdtraj remdtrajtemp 300.
        distance @10 @20
        '''
        input_file = 'data/Test_RemdTraj/traj.in'

        state_from_file = pt.load_cpptraj_file(input_file)
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


if __name__ == "__main__":
    unittest.main()
