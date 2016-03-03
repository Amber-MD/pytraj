import pytraj as pt

pdb = pt.load_pdb_rcsb("1l2y")
print(pdb)

h = pt.search_hbonds(pdb)
print(h)
print(h.donor_acceptor)
print('total solute hbonds: ', h.data['total_solute_hbonds'])
