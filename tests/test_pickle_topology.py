#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.compat import zip

def assert_equal_topology(top, new_top, traj):
    assert new_top.n_atoms == top.n_atoms, 'same n_atoms'
    assert new_top.n_residues == top.n_residues, 'same n_residues'
    assert new_top.n_mols == top.n_mols, 'same n_mols'
    # there are inverted bond indices [5292 5291] vs [5291 5292]
    # so use distance to assert
    aa_eq(pt.distance(traj, new_top.bond_indices), pt.distance(traj,
        top.bond_indices))
    # same for dihedral_indices
    aa_eq(pt.dihedral(traj, new_top.dihedral_indices), pt.dihedral(traj,
        top.dihedral_indices))
    aa_eq(new_top.dihedral_indices, top.dihedral_indices)
    aa_eq(new_top.mass, top.mass)
    aa_eq(new_top.charge, top.charge)

    assert [res.name for res in top.residues] == [res.name for res in
            new_top.residues], 'equal resnames'
    assert [atom.name for atom in top.atoms] == [atom.name for atom in
            new_top.atoms], 'equal resnames'

    for res, res_new in zip(top.residues, new_top.residues):
        assert res.first_atom_idx == res_new.first_atom_idx, 'first atom'
        assert res.last_atom_idx == res_new.last_atom_idx, 'last atom'


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
