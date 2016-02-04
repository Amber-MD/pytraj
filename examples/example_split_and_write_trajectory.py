import pytraj as pt
from pytraj.tools import split_and_write_traj

# use iterload for memory saving
traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

# split to 4 chunks
split_and_write_traj(traj, n_chunks=4, root_name='./output/trajx')

import os
os.system("ls ./output/trajx*")
