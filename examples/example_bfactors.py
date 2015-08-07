import pytraj as pt

pdb = pt.load_pdb_rcsb("1l2y")
bf = pt.bfactors(pdb, '@CA')
print(bf)
