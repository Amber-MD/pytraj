import os
import numpy as np
import unittest
import pytraj as pt
from pytraj.Topology import Topology
from pytraj.AtomMask import AtomMask
from pytraj.core.FileName import FileName
from pytraj.base import *
from pytraj.decorators import no_test
from pytraj.core.NameType import NameType

TRAJ = Trajectory("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")


class TestTopology(unittest.TestCase):
    def test_empty_top(self):
        top = Topology()
        assert top.is_empty() == True
        filename = "./data/Tc5b.top"
        top = Topology(filename)
        assert top.is_empty() == False

    def test_1(self):
        datadir = "./data/"
        filename = "./data/Tc5b.top"

        top = Topology(filename)
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
        top = Topology("./data/Tc5b.top")
        arr0 = top.atom_indices("@CA")

    def test_len(self):
        traj = Trajectory("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        assert len(top) == top.n_atoms

    def testLoadFromParmEd(self):
        import parmed as pmd
        fname = './data/Tc5b.top'

        orig_top = pt.load_topology(fname)
        parm = pmd.load_file('./data/Tc5b.top')
        top = pt.load_topology(parm)
        assert top.n_atoms == orig_top.n_atoms


if __name__ == "__main__":
    unittest.main()
