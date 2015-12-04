#!/usr/bin/env python

from __future__ import print_function
import sys
import unittest
import numpy as np
import subprocess
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.version import version


class TestExtraCoverage(unittest.TestCase):

    def test_extra_coverage(self):
        '''all kind of tests that do not belong to anywhere else
        '''
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        # show_versions
        pt.show_versions()
        pt._verbose()
        pt._verbose(False)
        print(version)

        # info
        pt.info()
        pt.info('parallel')
        pt.misc.parallel_info('pmap')
        pt.misc.parallel_info('openmp')
        pt.misc.parallel_info(None)

        eq([2, 3], [2, 3])
        # raise if comparing NaN
        self.assertRaises(ValueError, lambda: aa_eq(np.nan, np.nan))

        dslist = pt.multidihedral(traj)
        string_ = str(dslist[0])


if __name__ == "__main__":
    unittest.main()
