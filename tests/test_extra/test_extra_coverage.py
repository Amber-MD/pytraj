#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from utils import fn
from pytraj.utils import eq
from pytraj.testing import aa_eq
from pytraj.version import version
from pytraj.utils.get_common_objects import get_reference
from pytraj.utils import misc
from pytraj import all_actions, c_action

import pytest


class TestExtraCoverage(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

    def test_datafiles(self):
        traj = pt.datafiles.load_remd_ala2()
        assert len(traj.filelist) == 4, 'should have 4 replica trajs'
        filenames = [fn.split('/')[-1] for fn in traj.filelist]
        assert ['rem.nc.000', 'rem.nc.001', 'rem.nc.002', 'rem.nc.003'] == filenames, \
        'must have 4 replica trajs (rem.nc*)'

    def test_extra_coverage(self):
        '''all kind of tests that do not belong to anywhere else
        '''
        traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

        # show_versions
        pt.show_versions()
        pt._verbose()
        pt._verbose(False)

        # info
        pt.info()
        pt.info('parallel')
        misc.parallel_info('pmap')
        misc.parallel_info('openmp')
        misc.parallel_info(None)

        eq([2, 3], [2, 3])

        dslist = pt.multidihedral(traj)
        str(dslist[0])

    def testget_common_objects(self):
        # raises
        # raise if try to index traj()
        with pytest.raises(TypeError):
            get_reference(self.traj(), 3)
        with pytest.raises(TypeError):
            get_reference(self.traj(), None)

        # specify wrong mask
        with pytest.raises(TypeError):
            pt.superpose(self.traj[:], 3)

    def test_all_actions(self):
        with pytest.raises(ValueError):
            all_actions._assert_mutable(self.traj)

        with pytest.raises(AssertionError):
            # passing an Action object
            all_actions.do_action(
                self.traj, command='', action_class=c_action.Action_Angle())

        with pytest.raises(AssertionError):
            # passing a wrong class
            all_actions.do_action(
                self.traj, command='', action_class=pt.Trajectory)


if __name__ == "__main__":
    unittest.main()
