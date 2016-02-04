import pytraj as pt

traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")
top = traj.top

# get indices of CA atoms
print(pt.select("@CA", top))

# see how many H atoms
print(pt.select("@H=", top))
