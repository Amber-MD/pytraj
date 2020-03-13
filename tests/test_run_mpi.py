#!/usr/bin/env python

from __future__ import print_function
import subprocess
from glob import glob
import pytest

testlist = glob('test_*mpi/test_*py')

try:
    import mpi4py
except ImportError:
    mpi4py = None


@pytest.mark.skipif(mpi4py is None,
                    reason='require mpi4py')
@pytest.mark.parametrize('pyfile', testlist)
@pytest.mark.xfail(reason="fail on travis ci")
def test_all_mpi_scripts(pyfile):
    subprocess.check_call(['mpirun', '-n', '1', '--oversubscribe', 'python', pyfile])
