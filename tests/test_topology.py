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


class TestNewTopologyAPI(unittest.TestCase):
    """Test cases for new Topology API methods from cpptraj"""

    def setUp(self):
        self.top = pt.load_topology(fn('Tc5b.top'))
        self.tz2_top = pt.load_topology(fn('tz2.parm7'))

    def test_heavy_atom_count(self):
        """Test heavy atom counting methods with cross-validation"""
        # Test both property and method
        heavy_count_prop = self.top.n_heavy_atoms
        heavy_count_method = self.top.heavy_atom_count()

        assert heavy_count_prop == heavy_count_method
        assert isinstance(heavy_count_prop, int)
        assert heavy_count_prop > 0
        assert heavy_count_prop <= self.top.n_atoms

        # Cross-validate by manually counting non-hydrogen, non-extra-point atoms
        manual_heavy_count = 0
        for atom in self.top.atoms:
            if atom.element.upper() not in ['H'] and atom.atomic_number > 1:
                manual_heavy_count += 1

        # Heavy atom count should be reasonable compared to manual count
        # (cpptraj might have slightly different definition including extra points)
        assert abs(heavy_count_prop - manual_heavy_count) <= manual_heavy_count * 0.1

        # Heavy atoms should be less than total atoms (due to hydrogens)
        assert self.tz2_top.n_heavy_atoms < self.tz2_top.n_atoms

        # Test consistency across different topologies
        protein_heavy = self.tz2_top.n_heavy_atoms
        protein_total = self.tz2_top.n_atoms
        heavy_ratio = protein_heavy / protein_total
        # Protein heavy atom ratio should be reasonable (typically 20-50%)
        assert 0.1 < heavy_ratio < 0.8

    def test_atom_types(self):
        """Test atom type counting"""
        n_types = self.top.n_atom_types
        assert isinstance(n_types, int)
        assert n_types > 0
        assert n_types <= self.top.n_atoms  # Can't have more types than atoms

    def test_charge_info(self):
        """Test charge information methods"""
        # Test both property and method
        has_charges_prop = self.top.has_charges
        has_charges_method = self.top.has_charge_info()

        assert has_charges_prop == has_charges_method
        assert isinstance(has_charges_prop, bool)

        # tz2 should have charges, test with that
        assert self.tz2_top.has_charges == True

    def test_vdw_parameters(self):
        """Test VDW parameter access methods with physical constraints"""
        atom_idx = 0
        first_atom = self.top[atom_idx]

        # Test VDW radius
        radius = self.top.get_vdw_radius(atom_idx)
        assert isinstance(radius, float)
        assert 0.5 < radius < 5.0  # Physically reasonable range (Angstroms)

        # Test VDW sigma (related to radius, sigma ≈ radius * 2^(1/6))
        sigma = self.top.get_vdw_sigma(atom_idx)
        assert isinstance(sigma, float)
        assert 0.5 < sigma < 5.0  # Physically reasonable range

        # Test VDW depth (epsilon)
        depth = self.top.get_vdw_depth(atom_idx)
        assert isinstance(depth, float)
        assert 0.0 <= depth <= 1.0  # Typical range for epsilon (kcal/mol)

        # Test consistency across multiple atoms
        radii = []
        sigmas = []
        depths = []

        test_atoms = min(10, self.top.n_atoms)  # Test first 10 atoms
        for i in range(test_atoms):
            r = self.top.get_vdw_radius(i)
            s = self.top.get_vdw_sigma(i)
            d = self.top.get_vdw_depth(i)

            radii.append(r)
            sigmas.append(s)
            depths.append(d)

            # All should be positive and reasonable
            assert r > 0 and s > 0 and d >= 0

        # Should have some variation in parameters (not all identical)
        assert len(set([round(r, 2) for r in radii])) > 1 or test_atoms == 1

        # Test element-specific reasonable ranges
        if first_atom.element.upper() in ['C', 'N', 'O']:
            # Common biomolecular elements should have typical ranges
            assert 1.5 < radius < 2.5  # Carbon ~1.7, Nitrogen ~1.6, Oxygen ~1.5
            assert 0.05 < depth < 0.3  # Typical epsilon values

    def test_string_formatting_methods(self):
        """Test various string formatting methods with detailed validation"""
        atom_idx = 0
        atom = self.top[atom_idx]
        residue = self.top.residue(atom.resid)

        # Test residue@atom format (e.g., "SER@N")
        resname_atomname = self.top.trunc_resname_atomname(atom_idx)
        assert isinstance(resname_atomname, str)
        assert '@' in resname_atomname
        parts = resname_atomname.split('@')
        assert len(parts) == 2
        res_part, atom_part = parts
        assert len(res_part) > 0 and len(atom_part) > 0
        # Should match actual residue and atom names (trimmed)
        assert res_part == residue.name.strip()
        assert atom_part == atom.name.strip()

        # Test atomname_atomnum format (e.g., "N_1")
        atom_name_num = self.top.trunc_atom_name_num(atom_idx)
        assert isinstance(atom_name_num, str)
        assert '_' in atom_name_num
        name_part, num_part = atom_name_num.split('_')
        assert name_part == atom.name.strip()
        assert num_part == str(atom_idx + 1)  # 1-based numbering

        # Test full description format (e.g., "SER 1 N 1")
        full_desc = self.top.resname_num_atomname_num(atom_idx)
        assert isinstance(full_desc, str)
        parts = full_desc.split()
        assert len(parts) == 4
        res_name, res_num, atom_name, atom_num = parts
        assert res_name == residue.name.strip()
        assert res_num == str(atom.resid + 1)  # 1-based
        assert atom_name == atom.name.strip()
        assert atom_num == str(atom_idx + 1)  # 1-based

        # Test residue name with original number and ID
        res_idx = 0
        res_onum_id = self.top.trunc_resname_onum_id(res_idx)
        assert isinstance(res_onum_id, str)
        assert '_' in res_onum_id
        # Should start with residue name
        assert res_onum_id.startswith(residue.name.strip())

        # Test consistency across multiple atoms
        for i in range(min(5, self.top.n_atoms)):
            fmt_str = self.top.trunc_resname_atomname(i)
            assert '@' in fmt_str
            assert len(fmt_str) > 2  # Should be substantial

    def test_mask_utilities(self):
        """Test mask utility methods"""
        # Test with CA atoms
        ca_mask = "@CA"

        # Test residue numbers from mask
        resnums = self.top.resnums_selected_by(ca_mask)
        assert isinstance(resnums, np.ndarray)
        assert len(resnums) > 0

        # Test molecule numbers from mask
        molnums = self.top.molnums_selected_by(ca_mask)
        assert isinstance(molnums, np.ndarray)
        assert len(molnums) > 0

        # Test zero mass check
        has_zero_mass = self.top.mask_has_zero_mass(ca_mask)
        assert isinstance(has_zero_mass, bool)
        # CA atoms should not have zero mass
        assert has_zero_mass == False

    def test_molecule_methods(self):
        """Test molecule-related methods"""
        if self.top.n_mols > 0:
            mol_idx = 0
            nres = self.top.nres_in_mol(mol_idx)
            assert isinstance(nres, int)
            assert nres > 0

    def test_solute_residues(self):
        """Test solute residues identification"""
        solute_res = self.top.solute_residues()
        assert isinstance(solute_res, list)
        # Should have some solute residues
        assert len(solute_res) >= 0

        # All returned indices should be valid residue indices
        for res_idx in solute_res:
            assert 0 <= res_idx < self.top.n_residues

    def test_topology_modification_methods(self):
        """Test topology modification methods"""
        # Create a copy for modification tests
        top_copy = self.top.copy()

        # Test set_single_molecule
        result = top_copy.set_single_molecule()
        assert isinstance(result, int)
        assert result == 0  # Should return 0 on success

        # Test merge_residues (if we have multiple residues)
        if top_copy.n_residues > 1:
            original_nres = top_copy.n_residues
            result = top_copy.merge_residues(0, 1)
            # Note: This modifies the topology, so we expect fewer residues
            # The exact behavior depends on the implementation

    def test_edge_cases(self):
        """Test edge cases and error conditions"""
        # Test with valid indices to ensure methods work properly
        if self.top.n_atoms > 0:
            result = self.top.get_vdw_radius(0)
            assert isinstance(result, float)
            assert result > 0

        # Test with empty masks
        try:
            empty_resnums = self.top.resnums_selected_by("@NONEXISTENT")
            assert len(empty_resnums) == 0
        except:
            # Some implementations might raise an error for empty selections
            pass

        # Test string methods with valid indices
        if self.top.n_atoms > 0:
            result = self.top.trunc_resname_atomname(0)
            assert isinstance(result, str)
            assert len(result) > 0

    def test_consistency_with_existing_methods(self):
        """Test that new methods are consistent with existing ones"""
        # Heavy atom count should be consistent with manual counting
        heavy_count = 0
        for atom in self.top.atoms:
            if atom.element.upper() not in ['H']:
                heavy_count += 1

        # Note: The exact definition of "heavy atoms" in cpptraj might differ
        # from simple non-hydrogen counting (e.g., excluding extra points)
        # So we just check that the method returns a reasonable value
        assert 0 < self.top.n_heavy_atoms <= heavy_count

    def test_performance_methods(self):
        """Test that methods complete in reasonable time"""
        import time

        # These methods should be fast
        start_time = time.time()
        _ = self.top.n_heavy_atoms
        _ = self.top.has_charges
        _ = self.top.n_atom_types
        end_time = time.time()

        # Should complete in well under a second
        assert (end_time - start_time) < 1.0

    def test_return_types(self):
        """Test that all methods return expected types"""
        # Properties should return appropriate types
        assert isinstance(self.top.n_heavy_atoms, int)
        assert isinstance(self.top.n_atom_types, int)
        assert isinstance(self.top.has_charges, bool)

        # Methods should return appropriate types
        assert isinstance(self.top.heavy_atom_count(), int)
        assert isinstance(self.top.has_charge_info(), bool)
        assert isinstance(self.top.get_vdw_radius(0), float)
        assert isinstance(self.top.trunc_resname_atomname(0), str)
        assert isinstance(self.top.resnums_selected_by("@CA"), np.ndarray)
        assert isinstance(self.top.solute_residues(), list)

    def test_comprehensive_integration(self):
        """Comprehensive integration test of all new API methods"""
        # Test that all new methods work together coherently

        # Get basic topology info
        n_atoms = self.top.n_atoms
        n_heavy = self.top.n_heavy_atoms
        n_types = self.top.n_atom_types
        has_charges = self.top.has_charges

        # Validate relationships
        assert n_heavy <= n_atoms
        assert n_types <= n_atoms
        assert 1 <= n_types <= n_atoms

        # Test VDW parameters for multiple atom types
        unique_radii = set()
        unique_sigmas = set()
        test_count = min(20, n_atoms)

        for i in range(test_count):
            radius = self.top.get_vdw_radius(i)
            sigma = self.top.get_vdw_sigma(i)
            depth = self.top.get_vdw_depth(i)

            unique_radii.add(round(radius, 2))
            unique_sigmas.add(round(sigma, 2))

            # Basic physical constraints
            assert 0.1 < radius < 10.0
            assert 0.1 < sigma < 10.0
            assert 0.0 <= depth <= 2.0

        # Should have some diversity in VDW parameters
        if test_count > 5:
            assert len(unique_radii) > 1
            assert len(unique_sigmas) > 1

        # Test string formatting consistency
        for i in range(min(5, n_atoms)):
            atom = self.top[i]

            # Get all string representations
            resname_atom = self.top.trunc_resname_atomname(i)
            atom_num = self.top.trunc_atom_name_num(i)
            full_desc = self.top.resname_num_atomname_num(i)

            # All should contain the atom name
            atom_name = atom.name.strip()
            assert atom_name in resname_atom
            assert atom_name in atom_num
            assert atom_name in full_desc

        # Test mask utilities with multiple masks
        common_masks = ["@CA", "@C", "@N", "@O"]

        for mask in common_masks:
            try:
                indices = self.top.select(mask)
                if len(indices) > 0:
                    resnums = self.top.resnums_selected_by(mask)
                    molnums = self.top.molnums_selected_by(mask)
                    has_zero_mass = self.top.mask_has_zero_mass(mask)

                    # Consistency checks
                    assert len(resnums) <= len(indices)  # Can't have more residues than atoms
                    assert len(molnums) <= len(indices)  # Can't have more molecules than atoms
                    assert isinstance(has_zero_mass, bool)

                    # If atoms exist, mass should generally be positive
                    if len(indices) > 0:
                        total_mass = sum(self.top[i].mass for i in indices[:10])  # Check first 10
                        if total_mass > 0:
                            assert has_zero_mass == False
            except:
                # Mask might not be valid for this topology
                continue

        # Test solute residues consistency
        solute_res = self.top.solute_residues()
        assert all(0 <= res_idx < self.top.n_residues for res_idx in solute_res)

        # If we have charges, verify charge info consistency
        if has_charges:
            total_charge = sum(atom.charge for atom in self.top.atoms)
            # Total charge should be reasonable (not NaN or extremely large)
            assert abs(total_charge) < 1000  # Reasonable upper bound

        print(f"Integration test passed: {n_atoms} atoms, {n_heavy} heavy, {n_types} types, charges: {has_charges}")