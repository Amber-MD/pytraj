#!/usr/bin/env python

from __future__ import print_function
import sys
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import subprocess
from glob import glob


@unittest.skipIf(sys.platform == 'darwin', 'no test on osx for mpi')
class TestRunMPI(unittest.TestCase):

    def test_all_mpi_scripts(self):
        testlist = glob('test_mpi/test_*py')
        for testfile in testlist:
            print(testfile)
            subprocess.check_call(['mpirun', '-n', '4', 'python', testfile])


if __name__ == "__main__":
    unittest.main()
    # nosetests --with-coverage --cover-package pytraj -vs .
    # nosetests -vs --processes 6 --process-timeout 200 .
    # nosetests -vs --processes 6 --process-timeout 200 --with-coverage --cover-package pytraj .
