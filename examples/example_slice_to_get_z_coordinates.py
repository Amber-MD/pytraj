import pytraj as pt

# use iterload for memory saving
traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

# get z coordinates of CA atoms
# Note: very slow for large trajectory( GB or TB)

print(traj['@CA'].xyz[:, :, 2])
