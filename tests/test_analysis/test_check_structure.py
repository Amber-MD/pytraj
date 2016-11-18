#!/usr/bin/env python

from __future__ import print_function
import pytraj as pt
from utils import fn
from pytraj.testing import cpptraj_test_dir

tz2_fn = cpptraj_test_dir + '/tz2.parm7'
bad_pdb = cpptraj_test_dir + '/Test_CheckStructure/tz2.stretched.pdb'


def test_check_structure():
    traj = pt.iterload(bad_pdb, top=tz2_fn)
    command = 'offset 0.7'

    n_issues, out = pt.check_structure(traj, command)
    assert n_issues == 4
