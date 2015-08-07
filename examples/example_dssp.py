import pytraj as pt

pdb = pt.load_pdb_rcsb("1l2y")

out = pt.dssp(pdb)
print(out)
