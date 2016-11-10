#!/usr/bin/env python

from __future__ import print_function
import os
import unittest

import pytraj as pt
from pytraj.utils import aa_eq


leapin = """
source leaprc.protein.ff14SB
foo = sequence { ACE ALA NME }
saveamberparm foo foo.arm7 foo.crd
quit
"""

@unittest.skipIf(os.getenv("AMBERHOME") is None, "must set AMBERHOME")
class TestLeapNab(unittest.TestCase):
    tleap = os.getenv("AMBERHOME") + '/bin/tleap'
    antechamber = os.getenv("AMBERHOME") + '/bin/antechamber'

    @unittest.skipUnless(os.path.exists(tleap), 'no tleap found')
    def test_leap(self):
        traj = pt.io.load_leap(leapin)
        assert traj.n_atoms == 22
        assert traj.top.n_residues == 3

    @unittest.skipUnless(os.path.exists(antechamber), 'no antechamber found')
    def test_antechamber(self):
        fn = 'data/cro.ac'
        traj = pt.io.load_antechamber(fn)
        assert traj.n_atoms == 40
        assert traj.top.n_residues == 1
        aa_eq([atom.charge for atom in traj.top.atoms][:3], [-0.895800, 0.113200, 0.11110])
        aa_eq([atom.charge for atom in traj.top.atoms][-3:], [0.146500, 0.146500, 0.42300])


if __name__ == "__main__":
    unittest.main()
