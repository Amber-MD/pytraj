#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestAtomicCorr(unittest.TestCase):

    def test_atomiccorr(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        state = pt.load_batch(traj, '''
        atomiccorr out test.dat :1-13 byres
        ''')
        state.run()

        data = pt.atomiccorr(traj, ':1-13', byres=True)
        aa_eq(data, state.data[1].values)


if __name__ == "__main__":
    unittest.main()
