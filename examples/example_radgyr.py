import pytraj as pt

# use iterload for memory saving
traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

#
print(pt.radgyr(traj, '@CA', nomax=True))
