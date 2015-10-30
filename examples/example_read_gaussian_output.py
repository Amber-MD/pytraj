# Aim: read Gaussian output and write to pdb file as a Trajectory
# to view in VMD
# require: cclib
# install: "pip install cclib"

import pytraj as pt

traj = pt.tools.read_gaussian_output(filename="../tests/data/gaussian/GF2.log",
                                     top="../tests/data/gaussian/GF2.pdb")

# use options='model' to write multiple pdbs to a single file
pt.write_traj("output/traj_vmd.pdb", traj, options='model', overwrite=True)

# to view by VMD, just
# vmd ./output/traj_vmd.pdb
# you will see a movie with 56 frames.
