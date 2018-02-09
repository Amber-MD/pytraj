import pytraj as pt
import numpy as np

# turn on this to see cpptraj's stdout
# pt.set_cpptraj_verbose(True)

# use iterload for memory saving
traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

# for example you're only interested in a single Frame
# let's load it to memory
frame = traj[1]
print(frame)

# now you want to rotate omega dihedral angle for residue 3?

# make a list to store new frames
traj_xyz = []

for deg in range(-180, 180, 5):
    # since Frame does not hold Topology information, you need to
    # pass it
    pt.rotate_dihedral(frame, '3:omega:' + str(deg), top=traj.top)
    # always make a copy since
    # the `rotate_dihedral` inplace-change coords of frame
    traj_xyz.append(frame.xyz.copy())

# write multiple pdbs to a single file
# just like multiple pdbs when you load from rcsb
# you can use this pdb file to view in VMD without loading
# prmtop/gro/psf/...

new_traj = pt.Trajectory(xyz=np.asarray(traj_xyz), top=traj.top)
pt.write_traj(
    "./output/test0.pdb",
    new_traj,
    options='model',
    overwrite=True)

# you can use VMD to open the new file
# vmd ./output/test0.pdb
