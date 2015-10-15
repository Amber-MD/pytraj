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
        # Residue('ALA', resnum=0, icode=0, chainID=0)
        new_top = pt.Topology()
        print(top)

        for idx, (aname, atype, charge, mass, resnum, resname) in enumerate(zip(d['atom_name'],
                d['atom_type'], d['atom_charge'], d['atom_mass'], d['resnum'],
                d['resname'])):
            print(idx)
            atom = pt.core.Atom(aname, atype, charge, mass, resnum)
            residue = pt.core.Residue(resname, resnum)
            print(atom, residue)
            new_top.add_atom(atom, residue)
        print(new_top)


if __name__ == "__main__":
    unittest.main()
