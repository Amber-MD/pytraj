import pytraj as pt

# use `iterload` to save memory
traj = pt.iterload("../tests/data/tz2.ortho.nc",
                   "../tests/data/tz2.ortho.parm7")
print(traj)

# perform distance calculation with cpptraj mask
pt.distance(traj, ':2 :3')  # index starts from 1
pt.distance(traj, [':2 :3', '@CA @CB'])

# perform distance calculation with atom_indices array
pt.distance(traj, [[2, 5], [7, 9], [10, 13]])  # index starts from 0
