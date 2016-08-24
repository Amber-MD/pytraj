#!/usr/bin/env python
from __future__ import print_function

import sys
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj import utils
import doctest
from pytraj.externals.six import PY3
from pytraj import testing
from pytraj.datafiles import load_samples
from pytraj import energy_analysis
from pytraj.trajectory import frame
from pytraj import datafiles, cluster, nucleic_acid_analysis
from pytraj.analysis.c_action import actionlist
from pytraj.utils import convert
from pytraj.trajectory import frameiter
from pytraj import vector
from pytraj.datasets import datasetlist
from pytraj.analysis import base_holder
from pytraj.trajectory import trajectory_iterator
from pytraj.utils import check_and_assert
from pytraj.utils import get_common_objects
from pytraj.utils import decorators
from pytraj.analysis import nmr
from pytraj.datasets import array


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

    # @unittest.skipIf(sys.platform == 'darwin', 'linux testing only')
    def test_doc(self):

        modules = [convert, frameiter, vector, trajectory_iterator, ]
        if PY3:
            # avoid adding 'u' to string in PY2: u'GLU5_O-LYS8_N-H'
            if has_sander:
                modules.append(energy_analysis)
            additional_list = [
                frame,
                actionlist,
                cluster,
                datafiles,
                pt.all_actions,
                pt.topology,
                get_common_objects,
                pt.parallel.multiprocess,
                pt,
                pt.io,
                nucleic_acid_analysis,
                load_samples,
                pt.trajectory,
                decorators,
                pt.dssp_analysis,
                datasetlist,
                array,
                nmr,
                check_and_assert,
                pt.hbond_analysis,
                pt.tools,
                testing,
                utils,
                pt.matrix,
                base_holder,
            ]
            modules.extend(additional_list)
        assert get_total_errors(
            modules) == 0, 'doctest: failed_count must be 0'


if __name__ == "__main__":
    unittest.main()
