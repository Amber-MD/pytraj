import pytraj as pt

# use iterload for memory saving
traj = pt.iterload("../tests/data/md1_prod.Tc5b.x", "../tests/data/Tc5b.top")

# split to 4 chunks
traj.split_and_write_traj(n_chunks=4, root_name='./output/trajx')

import os
os.system("ls ./output/trajx*")
