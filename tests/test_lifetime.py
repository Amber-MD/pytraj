#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestLifeTime(unittest.TestCase):
    def test_lifetime_hbond_series(self):
        traj = pt.iterload("data/DPDP.nc", "data/DPDP.parm7")
        state = pt.load_pipeline(traj, '''
        hbond HB @N,H,C,O series
        # ’run’ is used here to process the trajectory and generate hbond data
        run
        # Perform lifetime analysis
        runanalysis lifetime HB[solutehb] out output/test.dat
        ''')
        state.run()
        print(state.data)


if __name__ == "__main__":
    unittest.main()
