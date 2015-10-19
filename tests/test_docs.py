#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import doctest
from pytraj.compat import PY3
from pytraj import testing


class TestDoc(unittest.TestCase):
    '''testing for light modules
    '''
    def test_doc(self):
        from pytraj.utils import convert
        def get_total_errors(modules):
            return sum([doctest.testmod(mod).failed for mod in modules])

        modules = [pt._get_common_objects,
                   pt._nastruct,
                   convert,
                  ]
        if PY3:
            # avoid adding 'u' to string in PY2: u'GLU5_O-LYS8_N-H'
            modules.append(pt.hbonds)
            # different formats between py2 and 3
            modules.append(pt.tools)
            modules.append(pt.parallel_mapping)
            modules.append(testing)
        assert get_total_errors(modules) == 0, 'doctest: failed_count must be 0'

if __name__ == "__main__":
    unittest.main()
