import pytraj as pt

traj = pt.iterload("../tests/data/md1_prod.Tc5b.x", "../tests/data/Tc5b.top")
top = traj.top

# get indices of CA atoms
print(top.select("@CA"))

# see how many H atoms
print(top.select("@H="))
