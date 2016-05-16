import pytraj as pt

# use load method to load all frames to memory
traj = pt.load("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

# compute center of geometry for residue 1, use all atoms
mask0 = ':1'
cog_0 = pt.center_of_geometry(traj, mask=mask0)
# should be ndarray, shape=(n_frames, 3) 
print(cog_0)

# compute center of geometry for residue 1, exclude H atoms
mask1 = ':1&!:1@H='
print(set(atom.name for atom in traj.top[mask1].atoms)) # should get {'OD1', 'CA', 'ND2', 'CB', 'O', 'C', 'N', 'CG'
cog_1 = pt.center_of_geometry(traj, mask=mask1)
print(cog_1)

# compute center of geometry for residue 1, 3, 5 and use all atoms
mask2 = ':1,3,5'
print(set(atom.name for atom in traj.top[mask2].atoms))
cog_2 = pt.center_of_geometry(traj, mask=mask2)
print(cog_2)

# compute center of geometry for residue 1, 3, 5 and use only CA atoms
mask3 = ':1,3,5@CA'
print(set(atom.name for atom in traj.top[mask3].atoms)) # should get {'@CA'}
cog_3 = pt.center_of_geometry(traj, mask=mask3)
print(cog_3)
