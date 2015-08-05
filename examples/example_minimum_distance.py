import pytraj as pt

pdb = pt.load_pdb_rcsb("1l2y")

# calculate minimum distance betwee two sets of atoms
out = pt.mindist(pdb, '@CA @H')
print (out)
