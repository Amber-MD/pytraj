#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn


class TestTopoloyObjects(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

    def test_atom(self):
        for idx, atom in enumerate(self.traj.top.atoms):
            assert atom.index == idx, 'residue index'

    def test_residue(self):
        for idx, res in enumerate(self.traj.top.residues):
            assert res.index == idx, 'residue index'


class TestSimplifiedTopology(unittest.TestCase):
    def test_simplified_topology(self):
        traj = pt.datafiles.load_trpcage()
        top = traj.top
        simplified_top = top.simplify()

        for res, simres in zip(top.residues, simplified_top.residues):
            assert res.name == simres.name
            assert res.index == simres.index

        for atom, simatom in zip(top.atoms, simplified_top.atoms):
            assert atom.name == simatom.name
            assert atom.atomic_number == simatom.atomic_number
            assert atom.charge == simatom.charge
            assert atom.type == simatom.type
            assert atom.mass == simatom.mass
            assert atom.element == simatom.element

        for simres in simplified_top.residues:
            for simatom in simres.atoms:
                assert simatom.resname == simres.name

            # test __iter__
            for simatom in simres:
                assert simatom.resname == simres.name


if __name__ == "__main__":
    unittest.main()
