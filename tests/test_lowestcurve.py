#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.c_analysis import c_analysis
from pytraj.datasets import CpptrajDatasetList

saved_data = '''
   0.348       9.3884
   0.548       9.2064
   0.748       9.1627
   0.948       9.3124
   1.148       9.3529
   1.348       9.3211
   1.548       9.2599
   1.748       9.5281
   1.948       9.4453
   2.148       9.5562
   2.348       9.6001
   2.548       9.7047
   2.748       9.5765
   2.948       9.4893
   3.148       9.6011
   3.348       9.9199
   3.548       9.8379
   3.748       9.8638
   3.948       9.9960
   4.148      10.1197
   4.348      10.2492
   4.548      10.2711
   4.748      10.3241
   4.948      10.0504
   5.148       9.9996
   5.348      10.2741
   5.548      10.4033
   5.748      10.6757
   5.948      10.7845
   6.148      10.9844
   6.348      11.3338
   6.548      12.1653
   6.748      13.1833
   6.948      13.4633
   7.148      13.2048
   7.348      13.8878
   7.548      13.6081
   7.748      14.8627
   7.948      14.0784
   8.148      14.7529
   8.348       0.0000
   8.548       0.0000
   8.748      13.2340'''

saved_data = np.array(
    [[x for x in line.split()] for line in saved_data.split('\n') if line],
    dtype='f4').T


class TestLowestCurve(unittest.TestCase):

    def test_lowestcurve_low_level(self):
        data = np.loadtxt('data/esurf_vs_rmsd.dat').T
        lc_data = pt.lowestcurve(data, points=10, step=0.2)
        aa_eq(saved_data, lc_data, decimal=3)


if __name__ == "__main__":
    unittest.main()
