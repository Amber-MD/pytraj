#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj import utils
import doctest
from pytraj.compat import PY3
from pytraj import testing

doctest.DONT_ACCEPT_BLANKLINE = False

def get_total_errors(modules):
    return sum([doctest.testmod(mod).failed for mod in modules])

class TestDoc(unittest.TestCase):
    '''testing for light modules
    '''
    def test_doc(self):
        from pytraj.utils import convert
        from pytraj import frameiter, vector, datasetlist
        from pytraj import trajectory_iterator
        from pytraj.parallel import pjob

        modules = [pt._get_common_objects,
                   pt._nastruct,
                   convert,
                   frameiter,
                   vector,
                   pjob,
                   datasetlist,
                   trajectory_iterator,
                  ]
        if PY3:
            # avoid adding 'u' to string in PY2: u'GLU5_O-LYS8_N-H'
            modules.append(pt.hbonds)
            # different formats between py2 and 3
            modules.append(pt.tools)
            modules.append(pt.parallel.parallel_mapping_multiprocessing)
            modules.append(testing)
            modules.append(utils)
        assert get_total_errors(modules) == 0, 'doctest: failed_count must be 0'

if __name__ == "__main__":
    unittest.main()
