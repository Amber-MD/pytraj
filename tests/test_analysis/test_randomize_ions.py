#!/usr/bin/env python

from __future__ import print_function
import os
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.testing import cpptraj_test_dir


class TestRandomizeIons(unittest.TestCase):
    def test_randomize_ions(self):
        fn = os.path.join(cpptraj_test_dir, 'Test_RandomizeIons',
                          'adh206.tip3p.rst7.gz')
        tn = os.path.join(cpptraj_test_dir, 'Test_RandomizeIons',
                          'adh206.ff10.tip3p.parm7.gz')
        saved_traj_name = os.path.join(cpptraj_test_dir, 'Test_RandomizeIons',
                                       'around.overlap.rst7.save')
        # Set default RNG back to Marsaglia
        pt.set_default_rng(0)
        traj = pt.iterload(fn, tn)
        traj_mut = traj[:]
        saved_traj = pt.iterload(saved_traj_name, traj.top)

        pt.randomize_ions(
            traj_mut,
            mask='@Na+',
            around=':1-16',
            by=5.0,
            overlap=3.0,
            seed=1)
        aa_eq(traj_mut.xyz, saved_traj.xyz, decimal=2)


if __name__ == "__main__":
    unittest.main()
