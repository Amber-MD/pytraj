#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.sandbox import to_parmed

try:
    import parmed as pmd
except ImportError:
    pmd = None


@unittest.skipIf(pmd is None, "Must install ParmEd")
class TestParmEdConverter(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
        self.parm = pmd.load_file(self.traj.top.filename)

    def parm_eq(self, parm, parm2):
        for atom, atom2 in zip(parm.atoms, parm2.atoms):
            assert atom.name == atom2.name, 'equal name'
            assert atom.type == atom2.type, 'equal type'
            assert atom.mass == atom2.mass, 'equal mass'
            assert atom.atomic_number == atom2.atomic_number, \
                'equal atomic_number'
            assert atom.residue.name == atom2.residue.name, 'residue name'
            aa_eq(atom.charge, atom2.charge, decimal=4)

    def test_parmed_converter(self):
        self.parm_eq(self.parm, to_parmed(self.traj))
        self.parm_eq(self.parm, self.traj.top.to_parmed())


if __name__ == "__main__":
    unittest.main()
