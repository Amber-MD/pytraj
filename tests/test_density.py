#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestDensity(unittest.TestCase):

    def test_density(self):
        fn = "data/DOPC.rst7"
        tn = "data/DOPC.parm7"

        traj = pt.load("data/DOPC.rst7", "data/DOPC.parm7")

        delta = '0.25'
        masks = [":PC@P31", ":PC@N31", ":PC@C2", ":PC | :OL | :OL2"]

        keys_no_space = [''.join(m.split()) for m in masks]

        mask_str = ' '.join(['"' + m + '"' for m in masks])
        density_types = ['number', 'mass', 'charge', 'electron']

        state_data_dict = dict()
        for density_type in density_types:
            command = """
            parm {parm}
            trajin {trajin}

            center ":PC | :OL | :OL2" origin
            density {density_type} delta {delta} {mask}
            """.format(parm=tn, trajin=fn,
                       delta=delta, mask=mask_str,
                       density_type=density_type)

            state = pt.load_cpptraj_state(command)
            state.run()

            state_data_dict[density_type] = state.data[1:].to_dict()

        pt.center(traj, '":PC | :OL | :OL2" origin')
        density_dict = {}
        for density_type in density_types:
            denstiy_data = pt.density(traj, mask=masks, delta=delta, density_type=density_type)
            density_dict[density_type] = denstiy_data

        # compate to cpptraj
        for density_type in density_types:
            for key in keys_no_space:
                aa_eq(state_data_dict[density_type][key],
                      density_dict[density_type][key])

        # assert raise: wrong density_type
        def func():
            pt.density(traj, mask=':WAT', density_type='hello')

        self.assertRaises(AssertionError, func)

if __name__ == "__main__":
    unittest.main()
