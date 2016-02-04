import pytraj as pt

traj = pt.load("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

# etract first 10 frames and write to CHARMM format
subtraj = traj[:10]
pt.write_traj(filename="./output/subtraj_0_CHARMM.dcd",
                traj=subtraj,
                top=traj.top,
                overwrite=True)

# make sure we can load the traj,
# use AMBER top is fine
charmm_traj = pt.load("./output/subtraj_0_CHARMM.dcd",
                        "../tests/data/Tc5b.top")

# calculate rmsd between old and saved traj for 1st frame
print(charmm_traj[0].rmsd(subtraj[0]))
assert charmm_traj[0].rmsd(subtraj[0]) < 1E-6

# another way
subtraj.save("./output/subtraj_1_CHARMM.dcd", overwrite=True)
charmm_traj_1 = pt.load("./output/subtraj_1_CHARMM.dcd",
                          "../tests/data/Tc5b.top")
print(charmm_traj_1[0].rmsd(subtraj[0]))
assert charmm_traj_1[0].rmsd(subtraj[0]) < 1E-6
