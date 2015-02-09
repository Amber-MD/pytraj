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
traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

# write to new format
traj.write("./output/md_charmm.dcd", overwrite=True)

# smart indexing, 3D array (n_frames, n_atoms, 3) for CA 
print (traj['@CA'])
# smart indexing, 3D array (n_frames, n_atoms, 3) for CA for 2th to 5th frame
print (traj[2:5]['@CA'])
# fancy indexing like Python/numpy array
# 3D array (n_frames, n_atoms, 3)
print (traj[2:5, :])

# internally get the data from cpptraj
# calculate radgyr for CA atoms using `traj`
d0 = calculate('radgyr', "@CA", traj)
print (d0[:])

# support reading trajectory file that cpptraj not yet handled
# hdf5
traj = mdio.load_hd5f("./data/ala2.h5")
# many more in github.com/pytraj/pytraj
