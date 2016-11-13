#!/usr/bin/env python

from __future__ import print_function
import subprocess
from glob import glob
import pytest

testlist = glob('test_*mpi/test_*py')

@pytest.mark.parametrize('pyfile', testlist)
def test_all_mpi_scripts(pyfile):
    subprocess.check_call(['mpirun', '-n', '4', 'python', pyfile])
