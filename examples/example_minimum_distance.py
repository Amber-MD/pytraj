import pytraj as pt

pdb = pt.load_pdb_rcsb("1l2y")

# calculate minimum distance betwee two sets of atoms
out = pt.mindist(pdb, '@CA @H')
print(out)

# use 2d array-like for two atom groups
out = pt.mindist(pdb, [[0, 3, 7], [6, 9, 100, 200]])
print(out)
