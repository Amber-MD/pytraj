import pytraj as pt

# Calculate density along a coordinate
# ask questio in Ambermailing list to know more about
# this calculation

traj = pt.iterload("../tests/data/DOPC.rst7", "../tests/data/DOPC.parm7")
print(traj)

print(pt.density(traj, 'delta 0.2 x mass @H='))

# get help from cpptraj
pt.info("density")
