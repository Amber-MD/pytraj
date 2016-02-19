#!/usr/bin/env python

from __future__ import print_function
import os
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestCheckOverlap(unittest.TestCase):

    def test_check_overlap(self):
        '''overlap checking
        '''
        tz2_bad = os.path.join(cpptraj_test_dir, 'Test_CheckStructure', 'tz2.stretched.pdb')
        traj = pt.iterload(tz2_bad, "data/tz2.parm7")
        aa_eq(pt.check_overlap(traj, options='offset 0.7'), [4,])


if __name__ == "__main__":
    unittest.main()
