#!/usr/bin/env python

from __future__ import print_function
import os
import unittest
import pytraj as pt
from pytraj.testing import aa_eq


DO_GIST = os.getenv("DO_GIST", False)

@unittest.skipUnless(DO_GIST, 'only do gist test if setting DO_GIST=True')
class TestGist(unittest.TestCase):
    def test_gist(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7", frame_slice=(0, 10))
        traj.autoimage('origin')

        command = "doorder doeij refdens 0.033422885325 gridcntr 1.44 0.67 0.29 \
                     griddim 10 12 10 gridspacn 2.0"
        data = pt.all_actions.gist(traj,
                do_order=True,
                do_eij=True,
                reference_density=0.033422885325,
                grid_center=(1.44, 0.67, 0.29), grid_dim=(10, 12, 10),
                grid_spacing=2.0,
                dtype='cpptraj_dataset')

        state_command = """
        parm data/tz2.ortho.parm7
        trajin data/tz2.ortho.nc
        autoimage origin
        gist {command}
        """.format(command=command)
        state = pt.load_cpptraj_state(state_command)
        state.run()

        data_dict = data.to_dict()
        data_state_dict = state.data[1:].to_dict()

        for key, state_key  in zip(sorted(data_dict.keys()), sorted(data_state_dict.keys())):
            aa_eq(data_dict[key], data_state_dict[state_key])

if __name__ == "__main__":
    unittest.main()
