#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq

cm = '''
pucker p1-as :1@C1' :1@C2' :1@C3' :1@C4' :1@O4'
pucker p2-as :2@C1' :2@C2' :2@C3' :2@C4' :2@O4'
pucker p3-as :3@C1' :3@C2' :3@C3' :3@C4' :3@O4'
pucker p1-cp :1@C1' :1@C2' :1@C3' :1@C4' :1@O4' cremer
pucker p2-cp :2@C1' :2@C2' :2@C3' :2@C4' :2@O4' cremer
pucker p3-cp :3@C1' :3@C2' :3@C3' :3@C4' :3@O4' cremer
'''


class TestPucker(unittest.TestCase):
    '''TestPucker
    '''

    def test_pucker(self):
        traj = pt.iterload('data/Test_NAstruct/adh026.3.pdb')
        state = pt.load_cpptraj_state(cm, traj)
        state.run()

        data_altona = pt.pucker(traj, resrange=range(3))
        data_cremer = pt.pucker(traj, resrange=range(3), method='cremer')
        aa_eq(data_altona.values, state.data[1:4].values)
        aa_eq(data_cremer.values, state.data[4:].values)

        # resrange=None
        data_full_residues = pt.pucker(traj, resrange=None)
        aa_eq(data_full_residues[:3].values, data_altona[:3].values)


if __name__ == "__main__":
    unittest.main()
