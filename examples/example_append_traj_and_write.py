from itertools import chain
import pytraj as pt
import os

traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

# write new traj with frames from 0 to 2, and 5 to 8
# skip last frame to follow python's convention.

# you can use traj[0:3], traj[5:9] but this does not save memory
# traj(0,3) like `range` in python3 or `xrange` in python2

pt.write_traj('traj_append.nc',
              traj=traj,
              frame_indices=chain(range(0, 3), range(5, 9)),
              top=traj.top,
              overwrite=True, )

# you can load many files to a single traj and write specific frames too
# traj = pt.iterload(['remd.x.000', 'remd.x.001', 'remd.x.009'], your_topology_file)
