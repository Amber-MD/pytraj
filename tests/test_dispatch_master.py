#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj import sandbox as sb


class TestSuperDispatch(unittest.TestCase):

    def test_radgyr_dispatch(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")

        # mask is an array-like
        aa_eq(sb._toy_radgyr(traj, [0, 3]), pt.radgyr(traj, '@1,4'))

        # frame_indices
        aa_eq(
            sb._toy_radgyr(traj,
                           '@CA',
                           frame_indices=[0, 3]),
            pt.radgyr(traj, '@CA')[[0, 3]])

        # frame_indices, mask in kwd
        aa_eq(
            sb._toy_radgyr(traj,
                           mask='@CA',
                           frame_indices=[0, 3]),
            pt.radgyr(traj, '@CA')[[0, 3]])

        # frame iterator
        aa_eq(
            sb._toy_radgyr(
                traj(0, 3),
                mask='@CA'),
            pt.radgyr(traj, '@CA')[[0, 1, 2]])


if __name__ == "__main__":
    unittest.main()
    # command: nosetests --with-coverage --cover-package pytraj -vs .
    # nosetests -vs --processes 6 --process-timeout 200 .
    # nosetests -vs --processes 6 --process-timeout 200 --with-coverage --cover-package pytraj .
