#!/usr/bin/env python

from __future__ import print_function
import pytraj as pt
from utils import fn
from pytraj.testing import cpptraj_test_dir, aa_eq

cpptraj_test_dir + '/tz2.parm7'
bad_pdb = cpptraj_test_dir + '/Test_CheckStructure/tz2.stretched.pdb'

parm_name = fn('DPDP.parm7')
traj_name = fn('DPDP.nc')
cm = "parm " + parm_name + "\ntrajin " + traj_name + "\ncheckchirality\n"


def test_check_structure():
    traj = pt.iterload(fn('DPDP.nc'), fn('DPDP.parm7'))
    out_dict = pt.check_chirality(traj)

    state = pt.load_cpptraj_state(cm)
    state.run()
    cpp_out_dict = state.data[1:].to_dict()

    for ((key, value), (cpp_key, cpp_value)) in zip(out_dict.items(),
                                                    cpp_out_dict.items()):
        aa_eq(value, cpp_value)
