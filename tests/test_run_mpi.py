#!/usr/bin/env python

from __future__ import print_function
import subprocess
from glob import glob

def test_all_mpi_scripts():
    testlist = glob('test_mpi/test_*py')
    for testfile in testlist:
        print(testfile)
        subprocess.check_call(['mpirun', '-n', '4', 'python', testfile])
