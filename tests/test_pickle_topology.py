#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        top = pt.iterload("./data/tz2.nc", "./data/tz2.parm7").top
        d = top.to_dict()

        # H = Atom('H', 'H', '0.0', '1.0', resnum=0)
        # res = Residue('ALA', resnum=0, icode=0, chainID=0)
        new_top = pt.Topology()

        for idx, (aname, atype, charge, mass, resnum, resname, mol_number) in enumerate(zip(d['atom_name'],
                d['atom_type'], d['atom_charge'], d['atom_mass'], d['resnum'],
                d['resname'], d['mol_number'])):
            atom = pt.core.Atom(aname, atype, charge, mass, resnum)
            residue = pt.core.Residue(resname, resnum)
            new_top.add_atom(atom, residue)

        new_top.add_bonds(d['bond_index'])
        new_top.add_dihedrals(d['dihedral_index'])
        
        assert new_top.n_atoms == top.n_atoms, 'same n_atoms'
        assert new_top.n_residues == top.n_residues, 'same n_residues'
        aa_eq(new_top.bond_indices, top.bond_indices)
        aa_eq(new_top.dihedral_indices, top.dihedral_indices)
        aa_eq(new_top.mass, top.mass)
        aa_eq(new_top.charge, top.charge)


if __name__ == "__main__":
    unittest.main()
