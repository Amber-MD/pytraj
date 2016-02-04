import pytraj as pt
import os

# use iterload for memory saving
traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

pt.write_traj("./output/trj.dcd", traj, overwrite=True)
os.system("ls ./output/trj.dcd")
