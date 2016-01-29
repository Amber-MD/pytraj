#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestMolsurf(unittest.TestCase):

    def setUp(self):
        self.traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

    def test_molsurf(self):
        traj = self.traj

        text = '''
        parm {0}
        trajin {1}
        molsurf @CA
        molsurf @CA probe 1.2
        molsurf @CA probe 1.2 offset 0.3
        '''.format(traj.top.filename, traj.filename)

        state = pt.load_cpptraj_state(text)
        state.run()
        cpp_data = state.data[1:].values

        atom_indices = traj.top.select("@CA")

        for mask in [atom_indices, '@CA']:
            aa_eq(pt.molsurf(traj, mask), cpp_data[0])
            aa_eq(pt.molsurf(traj, mask, probe=1.2), cpp_data[1])
            aa_eq(pt.molsurf(traj, mask, probe=1.2, offset=0.3), cpp_data[2])

    def test_surf(self):
        traj = pt.iterload('data/DPDP.nc', 'data/DPDP.parm7')

        text = '''
        parm {0}
        trajin {1}
        surf :1-12
        '''.format(traj.top.filename, traj.filename)

        state = pt.load_cpptraj_state(text)
        state.run()
        cpp_data = state.data[1:].values

        atom_indices = traj.top.select(":1-12")
        aa_eq(pt.surf(traj, ':1-12'), cpp_data)
        aa_eq(pt.surf(traj, atom_indices), cpp_data)


if __name__ == "__main__":
    unittest.main()
