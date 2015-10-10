#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.analyses.CpptrajAnalyses import Analysis_LowestCurve
from pytraj.datasets import CpptrajDatasetList


class TestLowestCurve(unittest.TestCase):
    @unittest.skip('not sure how this works yet')
    def test_lowestcurve(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        data = pt.rmsd(traj, mask='@CA')

        act = Analysis_LowestCurve()
        pt._verbose()
        dslist = CpptrajDatasetList()
        dslist.add_new('double', 'myrmsd')
        dslist[0].data = data
        act('points 10 myrmsd', dslist=dslist)


if __name__ == "__main__":
    unittest.main()
