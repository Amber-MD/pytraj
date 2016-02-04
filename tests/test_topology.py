import os
import numpy as np
import unittest
import pytraj as pt
from pytraj.compat import zip
from pytraj import Topology
from pytraj.core.c_core import AtomMask
from pytraj.base import *

TRAJ = Trajectory("./data/Tc5b.x", "./data/Tc5b.top")


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
        #
        top.strip("!@CA")
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

    def test_residue(self):
        for idx, res in enumerate(TRAJ.top.residues):
            assert idx == res.index, 'res.index'

    def test_get_iter(self):
        top = pt.load_topology("./data/DOPC.parm7")
        s = [atom.name for atom in top[":PC@H*"]]
        atom0 = top[":PC@H*"][0]

        old_natoms = top.n_atoms
        self.assertRaises(ValueError, lambda: top.join(top))
        top.join(top.copy())
        assert top.n_atoms == 2 * old_natoms

    def test_select_mask(self):
        top = pt.load_topology("./data/Tc5b.top")
        arr0 = top.atom_indices("@CA")

    def test_len(self):
        traj = Trajectory("./data/Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        assert len(top) == top.n_atoms

    def test_raise_RuntimeError(self):
        self.assertRaises(RuntimeError, lambda: pt.load_topology('dummy'))

    def test_get_atom_view(self):
        traj = pt.datafiles.load_ala3()
        top = traj.top
        atom_8 = top.atom(8)
        # set dummy number
        assert atom_8.gb_radius == 1.2, 'origin gb radius is 1.2'
        atom_8.gb_radius = 10.
        assert atom_8.gb_radius == 10, 'gb radius must be 10'
        assert top[8].gb_radius == 10, 'gb radius for top[10] must be 10'

    def test_join(self):
        # need to create traj0, traj1 to keep lifetime of t0, t1
        traj0 = pt.load_sample_data()
        t0 = traj0.top
        traj1 = pt.load_sample_data('tz2')
        t1 = traj1.top

        # +
        t2 = t0 + t1  # mimic ParmEd
        assert t2.n_atoms == t0.n_atoms + t1.n_atoms

        t0 += t1
        assert t0.n_atoms == t2.n_atoms

        # *

    def test_basic(self):
        '''slicing, select'''
        top = pt.load_topology("./data/Tc5b.top")

        #
        assert isinstance(top[0], Atom)
        assert isinstance(top[:2], pt.Topology)
        assert isinstance(top[:1], pt.Topology)
        assert isinstance(top[range(10)], pt.Topology)
        assert isinstance(top[list(range(10))], pt.Topology)
        assert isinstance(top[np.array(range(10))], pt.Topology)
        assert top[0].name == top['@1'][0].name

        # mask, AtomMask, python array, list
        atm = top("@CA")
        indices = atm.indices
        for a1, a2, a3, a4 in zip(top['@CA'], top[atm], top[indices],
                                  top[list(indices)]):
            assert a1.name == a2.name == a3.name == a4.name == 'CA'

        # check len
        assert len(top[:]) == top.n_atoms
        assert len(top[:10]) == 10

    def test_simplifed_topology(self):
        '''simplify'''
        top = pt.load_topology("./data/Tc5b.top")
        sim_top = top.simplify()

        for atom, sim_atom in zip(top.atoms, sim_top.atoms):
            assert atom.resname == sim_atom.resname, 'equal resname'
            assert atom.name == sim_atom.name, 'equal resname'
            assert atom.type == sim_atom.type, 'equal resname'
            assert atom.charge == sim_atom.charge, 'equal resname'
            assert atom.mass == sim_atom.mass, 'equal resname'

if __name__ == "__main__":
    unittest.main()
