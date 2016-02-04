import pytraj as pt

traj = pt.iterload('../tests/data/Tc5b.x', '../tests/data/Tc5b.top')

mat = pt.pairwise_rmsd(traj, '@CA')
print(mat)
