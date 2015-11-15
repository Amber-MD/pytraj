#!/usr/bin/env python

from __future__ import print_function
import sys
import unittest
import subprocess
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestExtraCoverage(unittest.TestCase):
    def test_extra_coverage(self):
        '''all kind of tests that do not belong to anywhere else
        '''
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        # show_versions
        pt.show_versions()
        pt._verbose()
        pt._verbose(False)

        # info
        pt.info()

        try:
            pt.to_mdtraj(traj)
        except ImportError:
            pass


if __name__ == "__main__":
    unittest.main()
