#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj import utils
import doctest
from pytraj.compat import PY3
from pytraj import testing
from pytraj.datafiles import load_samples
from pytraj.externals import energy
from pytraj import frame, datafiles, cluster, nucleic_acid_
from pytraj.c_action import actionlist


try:
    import sander
    has_sander = True
except ImportError:
    has_sander = False

doctest.DONT_ACCEPT_BLANKLINE = False


def get_total_errors(modules):
    return sum([doctest.testmod(mod).failed for mod in modules])


class TestDoc(unittest.TestCase):
    '''testing for light modules
    '''

    def test_doc(self):
        from pytraj.parallel import multiprocessing_
        modules = [pt.all_actions, ]
        if PY3:
            assert get_total_errors(
                modules) == 0, 'doctest: failed_count must be 0'


if __name__ == "__main__":
    unittest.main()
