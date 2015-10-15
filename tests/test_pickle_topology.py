#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import assert_equal_topology
from pytraj.compat import zip


class TestBuildAndPickleTopology(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

    def test_convert_to_dict_and_rebuild(self):
        '''test_convert_to_dict_and_rebuild
        '''
        top = self.traj.top
        d = top.to_dict()

        # H = Atom('H', 'H', '0.0', '1.0', resnum=0)
        # res = Residue('ALA', resnum=0, icode=0, chainID=0)
        new_top = pt.Topology()

        MOLNUM = 0

        for idx, (aname, atype, charge, mass, resnum, resname, mol_number) in enumerate(zip(d['atom_name'],
                d['atom_type'], d['atom_charge'], d['atom_mass'], d['resnum'],
                d['resname'], d['mol_number'])):
            atom = pt.core.Atom(aname, atype, charge, mass, resnum)
            atom.set_mol(mol_number)
            residue = pt.core.Residue(resname, resnum)
            if idx == 0:
                new_top.start_new_mol()
            if mol_number > MOLNUM:
                new_top.start_new_mol()
                MOLNUM += 1
            new_top.add_atom(atom, residue)

        new_top.add_bonds(d['bond_index'])
        new_top.add_dihedrals(d['dihedral_index'])

        assert_equal_topology(top, new_top, self.traj)

    def test_picle(self):
        '''test_picle
        '''
        pt.to_pickle(self.traj.top, 'output/new_top.pk')
        new_top = pt.read_pickle('output/new_top.pk')
        assert_equal_topology(self.traj.top, new_top, self.traj)

    def test_to_and_from_dict(self):
        cls = self.traj.top.__class__
        top = self.traj.top
        assert_equal_topology(top, cls.from_dict(top.to_dict()), self.traj)
        

if __name__ == "__main__":
    unittest.main()
