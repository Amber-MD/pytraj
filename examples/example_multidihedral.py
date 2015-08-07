import pytraj as pt

# use `iterload` to save memory
traj = pt.iterload("../tests/data/tz2.ortho.nc",
                   "../tests/data/tz2.ortho.parm7")

# search all possible dihedral types supported by cpptraj
dset = pt.multidihedral(traj)
print(dset)

# residue 3 to 7, skip every 2
dset2 = pt.multidihedral(traj, resrange='3-7')
print(dset2)

# residue 3 to 7, skip every 2
# use `range360`
dset = pt.multidihedral(traj, resrange='3-7', range360=True)
print(dset)
