# (require: mpi4py, numpy)
# mpirun -n 4 python mpi_cal_molsurf_0.py

# always add those lines to your code
import numpy as np
from mpi4py import MPI
import pytraj as pt
from pytraj.parallel import map_mpi as pymap
from pytraj.testing import aa_eq

comm = MPI.COMM_WORLD
# end. you are free to update anything below here

# split remd.x.000 to N cores and do calc_surf in parallel
root_dir = "../../tests/data/"
traj_name = root_dir + "tz2.ortho.nc"
parm_name = root_dir + "tz2.ortho.parm7"

# load to TrajectoryIterator
traj = pt.iterload(traj_name, parm_name)
#print(traj)

# mapping different chunk of `traj` in N cores
# need to provide `comm`
# save `total_arr` to rank=0
# others: total_arr = None
total_arr = pymap(pt.calc_molsurf, traj, "!:WAT", top=traj.top)

if comm.rank != 0:
    assert total_arr is None

if comm.rank == 0:
    # skip final array since its shape might be different from the rest
    t0 = np.asarray(total_arr[:-1]).flatten()
    t1 = np.asarray(total_arr[-1]).flatten()
    t = np.append(t0, t1)
    print(t)
    #print('total array len: ', t.shape[0])

    # assert to serial values
    t2 = pt.calc_molsurf(traj, "!:WAT", dtype='ndarray')
    aa_eq(t2, t)
