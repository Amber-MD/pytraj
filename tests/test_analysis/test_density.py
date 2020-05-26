#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
import pytest


class TestDensity(unittest.TestCase):
    def test_density(self):
        rst7 = fn('DOPC.rst7')
        tn = fn('DOPC.parm7')

        traj = pt.load(fn('DOPC.rst7'), fn('DOPC.parm7'))

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
            """.format(
                parm=tn,
                trajin=rst7,
                delta=delta,
                mask=mask_str,
                density_type=density_type)

            state = pt.load_cpptraj_state(command)
            state.run()

            state_data_dict[density_type] = state.data[1:].to_dict()

        pt.center(traj, '":PC | :OL | :OL2" origin')
        density_dict = {}
        for density_type in density_types:
            density_data = pt.density(
                traj, mask=masks, delta=delta, density_type=density_type)
            density_dict[density_type] = density_data

        # compate to cpptraj
        for density_type in density_types:
            for key in keys_no_space:
                aa_eq(state_data_dict[density_type][key],
                      density_dict[density_type][key])

        # assert raise: wrong density_type
        def func():
            pt.density(traj, mask=':WAT', density_type='hello')

        with pytest.raises(AssertionError):
            func()

        # test 'z' value
        saved_z_values = np.linspace(-24.1250, 23.8750, 193)
        aa_eq(density_data['z'], saved_z_values)


if __name__ == "__main__":
    unittest.main()
