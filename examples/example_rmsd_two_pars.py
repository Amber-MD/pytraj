"""calculate RMSD for two proteins with different topologies"""

from pytraj import io as mdio

# create trajectory objects (having frames)
traj_0 = mdio.load("../tests/data/Tc5b.crd", "../tests/data/Tc5b.top")
traj_1 = mdio.load("../tests/data/ala3.dcd", "../tests/data/ala3.psf")

# get new traj objects with given mask
traj_0_new = traj_0[":8-10@CA :frame"]
traj_1_new = traj_1[":1-3@CA :frame"]

# calculate rmsd between 1st frame of traj_1_new and 1st frame of traj_0_new
# best fit, no mass
print(traj_1_new[0].rmsd(traj_0_new[0]))
