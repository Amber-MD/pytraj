import pytraj as pt

# use iterload for memory saving
traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

# calculate phi residue 3
dset = pt.calc_phi(traj, resrange='3-7')

# calcuate autocorrelation function for 1st dataset
af = pt.acorr(dset[0])
print(af)
