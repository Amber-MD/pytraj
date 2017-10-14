#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.trajectory.trajectory_iterator import sort_filename_by_number
from mock import patch


class TestTrajectoryIterator(unittest.TestCase):
    def test_sorting_filelist(self):
        orig_list = [
            'md10.nc', 'md11.nc', 'md12.nc', 'md4.nc', 'md5.nc', 'md100.nc',
            'md6.nc', 'md7.nc', 'md8.nc', 'md9.nc'
        ]
        expected = [
            'md4.nc',
            'md5.nc',
            'md6.nc',
            'md7.nc',
            'md8.nc',
            'md9.nc',
            'md10.nc',
            'md11.nc',
            'md12.nc',
            'md100.nc',
        ]
        assert expected == sort_filename_by_number(
            orig_list), 'two filename list must be equal'

    def test_comprehensive(self):
        traj = pt.iterload(
            fn('Test_RemdTraj/rem.nc.000'),
            fn('Test_RemdTraj/ala2.99sb.mbondi2.parm7'))
        # temperature
        aa_eq(traj.temperatures, [
            300., 630.5, 630.5, 630.5, 630.5, 630.5, 630.5, 630.5, 492.2, 384.3
        ])

        # iterframe (already in doctest), just throwing raise to increase coverage score
        self.assertRaises(ValueError, lambda: traj.iterframe(rmsfit='crazy'))

        # raise
        # memory error if load larger than 1GB for xyz
        traj = pt.datafiles.load_tz2_ortho()
        for _ in range(2):
            traj._load(traj.filelist)

        with patch('pytraj.trajectory.trajectory_iterator.TrajectoryIterator._estimated_GB') as \
                mem:
             def mock_get(*args, **kwargs):
                 return 10000000
             mem.__get__ = mock_get
             self.assertRaises(MemoryError, lambda: traj.xyz)

        # can not find filename
        self.assertRaises(ValueError, lambda: traj._load(None))
        # has filename but does not have Topology
        self.assertRaises(
            ValueError,
            lambda: pt.TrajectoryIterator(fn('Test_RemdTraj/rem.nc.000'), top=None))
        self.assertRaises(
            ValueError,
            lambda: pt.TrajectoryIterator(fn('Test_RemdTraj/rem.nc.000')))
        # empty Topology
        self.assertRaises(
            ValueError,
            lambda: pt.TrajectoryIterator(fn('Test_RemdTraj/rem.nc.000'), top=pt.Topology()))
        # weird Topology
        self.assertRaises(
            ValueError,
            lambda: pt.TrajectoryIterator(fn('Test_RemdTraj/rem.nc.000'), top=pt.Frame))


if __name__ == "__main__":
    unittest.main()
    # nosetests --with-coverage --cover-package pytraj -vs .
    # nosetests -vs --processes 6 --process-timeout 200 .
    # nosetests -vs --processes 6 --process-timeout 200 --with-coverage --cover-package pytraj .
