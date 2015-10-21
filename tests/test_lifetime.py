#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.common_actions import lifetime


class TestLifeTime(unittest.TestCase):
    # TODO: assertion
    # not really understand how lifetime give the output.
    def test_lifetime_hbond_series(self):
        traj = pt.iterload("data/DPDP.nc", "data/DPDP.parm7")
        state = pt.load_pipeline(traj, '''
        hbond HB @N,H,C,O series
        run
        runanalysis lifetime HB[solutehb] out output/test.dat
        ''')
        state.run()
        hbonds_data = [d.values for d in state.data[2:21]]
        #exptected_lifetime_data = [d.values for d in state.data[-19:]]
        lifetime_data = [d.values for d in lifetime(hbonds_data, dtype='dataset')]

        #for arr0, arr1 in zip(exptected_lifetime_data, lifetime_data):
        #    aa_eq(arr0, arr1)


if __name__ == "__main__":
    unittest.main()
