import pytraj as pt

# use `iterload` to save memory
traj = pt.iterload("../tests/data/tz2.ortho.nc", "../tests/data/tz2.ortho.parm7")
print (traj)

# perform rmsd calculation to first frame, all atoms
pt.rmsd(traj, 0)

# perform rmsd calculation to last frame, all atoms but H
pt.rmsd(traj, -1, mask='!@H=')

# perform rmsd calculation to last frame, CA atoms
pt.rmsd(traj, -1, '@CA')
