import pytraj as pt

# to update coordinates, we need to load all data to memory
# so we use `load` method instead of `iterload`

traj = pt.load("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")
print(traj)

# nofit rmsd before fitting
print(pt.rmsd(traj, ref=-1, mask='@CA nofit'))

# superpose to last frame with mask='@CA'
traj.superpose(ref=-1, mask='@CA')

# make sure we DID superpose by calculate rmsd
print(pt.rmsd(traj, ref=-1, mask='@CA nofit'))
