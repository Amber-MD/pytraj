#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import doctest


class TestDoc(unittest.TestCase):
    def test_doc(self):
        failed_count = sum([doctest.testmod(pt._get_common_objects).failed,
                            doctest.testmod(pt.hbonds).failed,
                           ])
        assert failed_count == 0, 'doctest: failed_count must be 0'

if __name__ == "__main__":
    unittest.main()
