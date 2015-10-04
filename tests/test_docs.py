#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import doctest


class TestDoc(unittest.TestCase):
    def test_doct(self):
        doctest.testmod(pt._get_common_objects)

if __name__ == "__main__":
    unittest.main()
