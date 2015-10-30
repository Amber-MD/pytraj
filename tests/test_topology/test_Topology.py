import os
import numpy as np
import unittest
import pytraj as pt
from pytraj import Topology
from pytraj.core.cpp_core import AtomMask
from pytraj.base import *

TRAJ = Trajectory("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")


class TestTopology(unittest.TestCase):
    def test_empty_top(self):
        top = Topology()
        assert top.is_empty() == True
        filename = "./data/Tc5b.top"
        top = pt.load_topology(filename)
        assert top.is_empty() == False

    def test_1(self):
        datadir = "./data/"
        filename = "./data/Tc5b.top"

        top = pt.load_topology(filename)
        #top2 = top.modify_state_by_mask(AtomMask("!@CA"))
        #
        top.strip_atoms("!@CA")
        assert top.n_atoms == 20

        for atom in top.atoms:
            pass

        for res in top.residues:
            pass

        for mol in top.mols:
            pass

        for idx, atom in enumerate(top.atoms):
            pass
        assert idx + 1 == top.n_atoms

    def test_select_mask(self):
        top = pt.load_topology("./data/Tc5b.top")
        arr0 = top.atom_indices("@CA")

    def test_len(self):
        traj = Trajectory("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        assert len(top) == top.n_atoms

    def test_raise_RuntimeError(self):
        self.assertRaises(RuntimeError, lambda: pt.load_topology('dummy'))

    def test_get_atom_view(self):
        top = pt.datafiles.load_ala3().top
        atom_8 = top._atom(8)
        # set dummy number
        assert atom_8.gb_radius == 1.2, 'origin gb radius is 1.2'
        atom_8.gb_radius = 10.
        assert atom_8.gb_radius == 10, 'gb radius must be 10'
        assert top[8].gb_radius == 10, 'gb radius for top[10] must be 10'


if __name__ == "__main__":
    unittest.main()
