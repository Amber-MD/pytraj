#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestTopoloyObjects(unittest.TestCase):

    def setUp(self):
        self.traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

    def test_atom(self):
        for idx, atom in enumerate(self.traj.top.atoms):
            assert atom.index == idx, 'residue index'

    def test_residue(self):
        for idx, res in enumerate(self.traj.top.residues):
            assert res.index == idx, 'residue index'


if __name__ == "__main__":
    unittest.main()
