#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import numpy as np


class TestPairDist(unittest.TestCase):

    def test_general(self):
        traj = pt.iterload("./data/tz2.crd", "./data/tz2.parm7")

        for (mask, delta) in [('*', 0.1), ('@CA', '0.2')]:
            data = pt.pairdist(traj,
                               delta=delta,
                               mask=mask,
                               dtype='ndarray')
            data0 = data[0].T
            data1 = data[1].T

            txt = '''
            parm ./data/tz2.parm7
            trajin ./data/tz2.crd
            pairdist out test.out mask {0} delta {1}
            '''.format(mask, delta)

            # take data, skip DatasetTopology
            cpp_data = pt.datafiles.load_cpptraj_output(txt)[1:].values
            cpp_distance, cpp_Pr = cpp_data[0].T
            _, cpp_std = cpp_data[1].T

            aa_eq(data0[0], cpp_distance)  # distance
            aa_eq(data1[0], cpp_distance)  # distance
            aa_eq(data0[1], cpp_Pr, decimal=2)  # Pr
            aa_eq(data1[1], cpp_std, decimal=2)  # std


if __name__ == "__main__":
    unittest.main()
