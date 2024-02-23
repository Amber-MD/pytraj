#!/usr/bin/env python

import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq

cm = '''
multipucker
'''

def test_multipucker():
    traj = pt.iterload(fn('Test_NAstruct/adh026.3.pdb'))
    state = pt.load_cpptraj_state(cm, traj)
    state.run()

    print(state.data[1:])

    data = pt.multipucker(traj)
    print(data)
