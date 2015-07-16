import pytraj as pt

pdb = pt.load_pdb_rcsb("1l2y")
print (pdb)

dslist = pt.search_hbonds(pdb)
print (dslist)
print (dslist.keys())
