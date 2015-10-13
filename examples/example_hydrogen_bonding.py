import pytraj as pt

pdb = pt.load_pdb_rcsb("1l2y")
print(pdb)

dslist = pt.search_hbonds(pdb)
print(dslist)
print(dslist.donor_aceptor)
print('total solute hbonds: ', dslist.data['total_solute_hbonds'])
