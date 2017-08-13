from __future__ import print_function
import sys
import pytraj as pt
import numpy as np

try:
    import sander
    from parmed.utils.netcdf import netcdf_file
except ImportError:
    print("Example requires pysander and  parmed installed")
    sys.exit(0)

fh = netcdf_file(
    'data/mdfrc', mmap=False)  # change mdfrc to your force filename
forces = fh.variables['forces']

# do similar thing for velocities

top = pt.load_topology(
    'data/mdfrc.prmtop')  # change prmtop to your parm7 filename
traj = pt.Trajectory(xyz=np.ascontiguousarray(forces[:], dtype='f8'), top=top)
print('original trajectory', traj)

# strip atom atoms
# change ':1-100' to your mask (e.g ':WAT')
new_traj = pt.strip(':1-100', traj)
print('stripped trajectory', new_traj)

# write trajectory
new_traj.save('test.nc', overwrite=True)
