import numpy as np
import pytraj as pt
import os

traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")
ref = pt.iterload("../tests/data/Tc5b.nat.crd", traj.top)

data = pt.rmsd(traj, ref=ref, mask='@CA')
print(data)

pt.write_traj("test_indices.pdb", traj,
              frame_indices=[2, 8, 9],
              overwrite=True)

indices = np.where(data < 5.0)[0]
# filter, write only frames having rmsd < 5.0
pt.write_traj("test_indices_2.pdb", traj,
              frame_indices=indices,
              overwrite=True)

# make sure we did save the traj
os.system("ls test_indices.pdb")
os.system("ls test_indices_2.pdb")
