#!/usr/bin/env python

from __future__ import print_function
import os
import unittest

import pytraj as pt
from pytraj.testing import aa_eq

leapin = """
source leaprc.protein.ff14SB
foo = sequence { ACE ALA NME }
saveamberparm foo foo.arm7 foo.crd
quit
"""

tleap = os.getenv("AMBERHOME", '') + '/bin/tleap'

@unittest.skipUnless(os.path.exists(tleap), 'no tleap found')
def test_leap():
    traj = pt.io.load_leap(leapin)
    assert traj.n_atoms == 22
    assert traj.top.n_residues == 3
