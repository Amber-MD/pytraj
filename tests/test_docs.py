#!/usr/bin/env python
from __future__ import print_function

import unittest
import pytraj as pt
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
from pytraj.trajectory import trajectory, trajectory_iterator
from pytraj.utils import check_and_assert
from pytraj.utils import get_common_objects
from pytraj.utils import decorators
from pytraj.analysis import nmr
from pytraj.datasets import array
from pytraj.builder import build


try:
    has_sander = True
except ImportError:
    has_sander = False

doctest.DONT_ACCEPT_BLANKLINE = False


def get_total_errors(modules):
    return sum([doctest.testmod(mod).failed for mod in modules])


@unittest.skipUnless(PY3, 'only test for PY3')
class TestDoc(unittest.TestCase):
    '''testing for light modules
    '''

    # @unittest.skipIf(sys.platform == 'darwin', 'linux testing only')
    def test_doc(self):

        modules = [convert, vector]
        # avoid adding 'u' to string in PY2: u'GLU5_O-LYS8_N-H'
        if has_sander:
            modules.append(energy_analysis)
        additional_list = [
            frame,
            actionlist,
            datafiles,
            pt.topology,
            get_common_objects,
            pt.parallel.multiprocess,
            pt,
            nucleic_acid_analysis,
            load_samples,
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
        assert not get_total_errors(modules)

    def test_doc_trajectory(self):
        modules = [trajectory, trajectory_iterator, frameiter]
        assert not get_total_errors(modules)

    def test_doc_io(self):
        modules = [pt.io]
        assert not get_total_errors(modules)

    def test_doc_all_actions(self):
        modules = [pt.all_actions,]
        assert not get_total_errors(modules)

    def test_builder(self):
        modules = [build,]
        assert not get_total_errors(modules)

    def test_clustering(self):
        modules = [cluster,]
        assert not get_total_errors(modules)

if __name__ == "__main__":
    unittest.main()
