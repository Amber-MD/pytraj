#!/usr/bin/env python

from __future__ import print_function
import os
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.utils.context import capture_stdout

# local
from utils import fn, tz2_ortho_trajin, tz2_ortho_top

DO_GIST = os.getenv("DO_GIST", False)


@unittest.skipUnless(DO_GIST, 'only do gist test if setting DO_GIST=True')
def test_gist():
    traj = pt.iterload(tz2_ortho_trajin, tz2_ortho_top, frame_slice=(0, 10))
    traj.autoimage('origin')

    command = "doorder doeij refdens 0.033422885325 gridcntr 1.44 0.67 0.29 \
                 griddim 10 12 10 gridspacn 2.0"

    data = pt.all_actions.gist(
        traj,
        do_order=True,
        do_eij=True,
        reference_density=0.033422885325,
        grid_center=(1.44, 0.67, 0.29),
        grid_dim=(10, 12, 10),
        grid_spacing=2.0,
        dtype='cpptraj_dataset')

    state_command = """
    parm {}
    trajin {}
    autoimage origin
    gist {}
    """.format(tz2_ortho_top, tz2_ortho_trajin, command)
    state = pt.load_cpptraj_state(state_command)
    with capture_stdout() as (out, _):
        state.run()

    data_dict = data.to_dict()
    data_state_dict = state.data[1:].to_dict()

    for key, state_key in zip(
            sorted(data_dict.keys()), sorted(data_state_dict.keys())):
        aa_eq(data_dict[key], data_state_dict[state_key])
