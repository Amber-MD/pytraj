#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/tz2.crd", "./data/tz2.parm7")
        txt = '''
        parm data/tz2.parm7
        trajin data/tz2.crd
        pairdist mask "*" delta 0.1
        '''
        cpp_data = pt.datafiles.load_cpptraj_output(txt)
        print(cpp_data)
        pt._verbose()
        print(pt.pairdist(traj, delta=0.1))


if __name__ == "__main__":
    unittest.main()
