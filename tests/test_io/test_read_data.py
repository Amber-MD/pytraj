#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.testing import cpptraj_test_dir


class TestReadWriteData(unittest.TestCase):
    def test_read(self):
        fn0 = "{cpptraj_test_dir}/Test_Vector/vtest.dat.6.save".format(
            cpptraj_test_dir=cpptraj_test_dir)
        fn1 = "{cpptraj_test_dir}/Test_Vector/vtest.dat.7.save".format(
            cpptraj_test_dir=cpptraj_test_dir)

        cm = """
        readdata {} vector name v6
        readdata {} vector name v7
        """.format(fn0, fn1)

        state = pt.load_cpptraj_state(cm)
        state.run()

        data0 = pt.io.read_data(fn0, options='vector name v6')
        data1 = pt.io.read_data(fn1, options='vector name v7')
        aa_eq(state.data[0].values, data0.values)
        aa_eq(state.data[1].values, data1.values)


if __name__ == "__main__":
    unittest.main()
