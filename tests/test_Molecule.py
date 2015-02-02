from pytraj.Molecule import Molecule

mol = Molecule(1, 304)

print(mol.begin_atom)
print(mol.end_atom)
print(mol.n_atoms)
print(mol.is_solvent())


# test wrong inputs
# mol2 = Molecule(1, 2, 3)
