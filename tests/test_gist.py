#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestGist(unittest.TestCase):

    def test_gist(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        command = "doorder doeij refdens 0.033422885325 gridcntr 1.44 0.67 0.29 griddim 10 12 10 gridspacn 2.0"

        state_command = """
        parm data/tz2.ortho.parm7
        trajin data/tz2.ortho.nc 1 10
        autoimage origin
        gist {}
        """.format(command)

        state = pt.load_cpptraj_state(state_command)
        state.run()
        print(state.data[1:].to_dict())

if __name__ == "__main__":
    unittest.main()
