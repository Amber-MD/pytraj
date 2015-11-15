#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestMultiVector(unittest.TestCase):

    def test_multivector(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        state = pt.load_batch(traj, '''
        multivector resrange 3-7 name1 C name2 N
        ''')
        state.run()
        cpp_data = state.data[1:].values

        aa_eq(
            pt.multivector(traj,
                           resrange='3-7',
                           names=('C', 'N'),
                           dtype='ndarray'),
            cpp_data)
        aa_eq(
            pt.multivector(traj,
                           resrange='3-7',
                           names='C N',
                           dtype='ndarray'),
            cpp_data)
        aa_eq(
            pt.multivector(traj,
                           resrange='3-7',
                           names='name1 C name2 N',
                           dtype='ndarray'),
            cpp_data)


if __name__ == "__main__":
    unittest.main()
