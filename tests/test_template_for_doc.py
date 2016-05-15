#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
import doctest
from pytraj.compat import PY3
from pytraj import testing


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
        modules = [pt.frame, ]
        if PY3:
            assert get_total_errors(
                modules) == 0, 'doctest: failed_count must be 0'


if __name__ == "__main__":
    unittest.main()
