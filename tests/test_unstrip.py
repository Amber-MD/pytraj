#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestUnstrip(unittest.TestCase):

    def test_unstrip(self):
        from pytraj.datasets import CpptrajDatasetList

        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        dslist = CpptrajDatasetList()
        actlist = pt.ActionList(
            ['strip !@CA', 'unstrip', 'radgyr nomax'],
            top=traj.top,
            dslist=dslist)

        actlist.compute(traj)

        # make sure that after stripping and unstrip the Frame coords are restored
        aa_eq(pt.radgyr(traj), dslist.values)

        dslist2 = CpptrajDatasetList()
        actlist2 = pt.ActionList(
            ['strip !@CA', 'radgyr nomax'],
            top=traj.top,
            dslist=dslist2)
        actlist2.compute(traj)
        # make sure get correct radgyr after stripping all but CA atoms
        aa_eq(dslist2.values, pt.radgyr(traj, '@CA'))


if __name__ == "__main__":
    unittest.main()
