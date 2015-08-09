# github.com/pytraj/pytraj
# pytraj: API for cpptraj: extending flexiblity for data analysis

# (install the most updated pytraj version on github)
# copy this file to ./tests folder in pytraj home
# cd tests/
# python ./for_Amber_meeting_2015.py

from pytraj import io as mdio
from pytraj import calculate

# smart loading
# topology = mdio.load("./data/Tc5b.top")
# if traj file exists, load to traj with topology
# can handle list of trajectory files too
traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

# write to new format
traj.write("./output/md_charmm.dcd", overwrite=True)

# smart indexing, 3D array (n_frames, n_atoms, 3) for CA
print(traj['@CA'])
# smart indexing, 3D array (n_frames, n_atoms, 3) for CA for 2th to 5th frame
print(traj[2:5]['@CA'])
# fancy indexing like Python/numpy array
# 3D array (n_frames, n_atoms, 3)
print(traj[2:5, :])

# expose cpptraj Frame object to numpy/scipy libray
import numpy as np
frame0 = traj[0]
# all calculation with arr0 will update frame0
arr0 = np.asarray(frame0.buffer2d[:])
print(arr0.shape)

# internally get the data from cpptraj
# calculate radgyr for CA atoms using `traj`
d0 = calculate('radgyr', traj, "@CA")
print(d0[:])

# calculate DSSP for first 3 frames, return either 'int' or 'str' array
from pytraj.common_actions import calc_dssp
arr0 = calc_dssp(traj[:3], '@CA', dtype='int')
arr1 = calc_dssp(traj[:3], '@CA', dtype='str')
print(arr0)
print(arr1)

# support reading trajectory file that cpptraj not yet handled
# hdf5
traj = mdio.load_hdf5("./data/ala2.h5")
# many more in github.com/pytraj/pytraj
