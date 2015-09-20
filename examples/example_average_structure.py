import pytraj as pt

pdb = pt.load_pdb_rcsb("1l2y")

avg_pdb = pt.mean_structure(pdb, '@CA')
print(avg_pdb)
