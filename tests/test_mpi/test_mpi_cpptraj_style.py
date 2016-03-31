# always add those lines to your code
import sys
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq

try:
    from mpi4py import MPI
except ImportError:
    print('skip if not having mpi4py')
    sys.exit(0)

comm = MPI.COMM_WORLD
# end. you are free to update anything below here

# split remd.x.000 to N cores and do calc_surf in parallel
root_dir = "data/"
traj_name = root_dir + "tz2.ortho.nc"
parm_name = root_dir + "tz2.ortho.parm7"

# load to TrajectoryIterator
traj = pt.iterload(traj_name, parm_name)

# save `total_arr` to rank=0
# others: total_arr = None
total_arr = pt.pmap_mpi(
    ['autoimage', 'center :2', 'distance :3 :7', 'angle :3 :7 :8'], traj)

if comm.rank != 0:
    assert total_arr is None

if comm.rank == 0:
    # assert to serial
    from pytraj.tools import dict_to_ndarray
    arr = dict_to_ndarray(total_arr)

    t0 = pt.center(traj[:].autoimage(), ':2')
    aa_eq(pt.distance(t0, ':3 :7'), arr[0])
    aa_eq(pt.angle(t0, ':3 :7 :8'), arr[1])
