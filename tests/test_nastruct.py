#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt

class TestNastruct(unittest.TestCase):
    def test_nupars(self):
        fn = "./data/Test_NAstruct/adh026.3.pdb"
        traj = pt.iterload(fn, fn)
        dslist = pt.nastruct(traj)

        text = '''
        parm "./data/Test_NAstruct/adh026.3.pdb"
        trajin "./data/Test_NAstruct/adh026.3.pdb"
        nastruct
        '''

        state = pt.load_cpptraj_state(text)
        state.run()
        print(state.data)


if __name__ == "__main__":
    unittest.main()
