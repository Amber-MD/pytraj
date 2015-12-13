#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.sandbox import to_parmed
import parmed as pmd


class TestParmEdConverter(unittest.TestCase):

    def test_parmed_converter(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        parm = pmd.load_file(traj.top.filename)
        parm2 = to_parmed(traj)

        for atom, atom2 in zip(parm.atoms, parm2.atoms):
            assert atom.name == atom2.name, 'equal name'
            assert atom.type == atom2.type, 'equal type'
            assert atom.mass == atom2.mass, 'equal mass'
            assert atom.atomic_number == atom2.atomic_number, 'equal atomic_number'
            assert atom.residue.name == atom2.residue.name, 'residue name'
            aa_eq(atom.charge, atom2.charge, decimal=4)

if __name__ == "__main__":
    unittest.main()
