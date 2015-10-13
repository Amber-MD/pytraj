#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import doctest


class TestDoc(unittest.TestCase):
    def test_doc(self):
        def get_total_errors(modules):
            return sum([doctest.testmod(mod).failed for mod in modules])

        modules = [pt._get_common_objects,
                   pt.hbonds,
                  ]
        assert get_total_errors(modules) == 0, 'doctest: failed_count must be 0'

if __name__ == "__main__":
    unittest.main()
