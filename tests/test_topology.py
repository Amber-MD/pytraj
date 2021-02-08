import numpy as np
import unittest
import pytraj as pt
from pytraj import Topology, Trajectory, Atom
from pytraj.core.elements import mass_atomic_number_dict
from pytraj.testing import aa_eq

# local
from utils import fn
import pytest

TRAJ = Trajectory(fn("Tc5b.x"), fn("Tc5b.top"))

tc5b_top = fn('Tc5b.top')


class TestTopology(unittest.TestCase):
    def test_empty_top(self):
        top = Topology()
        assert top.is_empty() == True
        top = pt.load_topology(tc5b_top)
        assert top.is_empty() == False

    def test_1(self):
        top = pt.load_topology(tc5b_top)
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
        top = pt.load_topology(fn('DOPC.parm7'))
        [atom.name for atom in top[":PC@H*"]]
        top[":PC@H*"][0]

        old_natoms = top.n_atoms
        with pytest.raises(ValueError):
            top.join(top)
        top.join(top.copy())
        assert top.n_atoms == 2 * old_natoms

    def test_select_mask(self):
        top = pt.load_topology(tc5b_top)
        top.atom_indices("@CA")

    def test_len(self):
        traj = Trajectory(fn("Tc5b.x"), tc5b_top)
        top = traj.top
        assert len(top) == top.n_atoms

    def test_raise_RuntimeError(self):
        with pytest.raises(RuntimeError):
            pt.load_topology('dummy')

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
        top = pt.load_topology(fn('Tc5b.top'))

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

        # API
        top.bond_indices
        top.angle_indices
        top.dihedral_indices

    def test_simplifed_topology(self):
        '''simplify'''
        top = pt.load_topology(fn('Tc5b.top'))
        sim_top = top.simplify()

        assert sim_top.select('@CA').tolist() == top.select('@CA').tolist()

        for atom, sim_atom in zip(top.atoms, sim_top.atoms):
            assert atom.resname == sim_atom.resname
            assert atom.name == sim_atom.name
            assert atom.type == sim_atom.type
            assert atom.charge == sim_atom.charge
            assert atom.mass == sim_atom.mass

        # API
        atom = sim_top.atoms[0]
        atom.residue
        atom.residue.name
        atom.residue.index
        atom.bond_partners


def test_mass_atomic_number_dict():
    top = pt.load_topology(fn("tz2.parm7"))
    mass_list = []

    for atom in top:
        mass_list.append(mass_atomic_number_dict[atom.atomic_number])
    aa_eq(mass_list, top.mass, decimal=2)


def test_toplogy_mass():
    traj = pt.iterload(fn("Tc5b.x"), fn("Tc5b.top"))
    top = traj.top

    mlist = []
    for atom in top:
        mlist.append(atom.mass)
    mlist = np.array(mlist)
    aa_eq(top.mass, mlist)
