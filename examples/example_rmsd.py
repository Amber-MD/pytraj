import pytraj as pt

# use `iterload` to save memory
traj = pt.iterload("../tests/data/tz2.ortho.nc",
                   "../tests/data/tz2.ortho.parm7")
print(traj)

# perform rmsd calculation to first frame, all atoms
pt.rmsd(traj, ref=0)

# perform rmsd calculation to last frame, all atoms but H
pt.rmsd(traj, mask='!@H=', ref=-1)

# perform rmsd calculation to last frame, CA atoms
pt.rmsd(traj, '@CA', ref=-1)

# perform rmsd calculation to stripped-atom reference
pt.rmsd(traj(mask='@CA'), ref=traj[2:3, '@CA'])
# `ref=traj[2:3, '@CA']` is equal to `ref=traj[2:3]['@CA'][0]`
