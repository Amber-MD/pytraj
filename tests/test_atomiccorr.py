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

    def test_issue_172_cpptraj(self):
        '''Internal Error: Attempting to assign to Frame with external memory
        '''
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        arr_out_of_memory = pt.atomiccorr(traj(0, 8, 2), '@CA')
        arr_in_memory = pt.atomiccorr(traj[0: 8: 2], '@CA')
        aa_eq(arr_out_of_memory, arr_in_memory)


if __name__ == "__main__":
    unittest.main()
