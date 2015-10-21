#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.analyses.CpptrajAnalyses import Analysis_LowestCurve
from pytraj.datasets import CpptrajDatasetList


class TestLowestCurve(unittest.TestCase):
    #@unittest.skip('not sure how this works yet')
    def test_lowestcurve(self):
        data = np.loadtxt('data/esurf_vs_rmsd.dat').T
        act = Analysis_LowestCurve()
        dslist = CpptrajDatasetList()

        # add rmsd
        dslist.add_new('double', 'myrmsd')
        dslist[0].data = data[0]

        ## add esurf
        dslist.add_new('double', 'mysurf')
        dslist[1].data = data[1]

        #dslist.add_new('xymesh', 'mysurf')
        #dslist[0]._append_from_array(data.T)

        #act('points 10 mysurf', dslist=dslist)
        act('points 10 myrmsd mysurf', dslist=dslist)
        print([d for d in dslist])


if __name__ == "__main__":
    unittest.main()
