# (require: mpi4py, numpy)
# mpirun -n 4 python mpi_cal_molsurf_0.py

# always add those lines to your code
import numpy as np
from mpi4py import MPI
from pytraj.parallel import map as pymap
from pytraj import io
import pytraj.common_actions as pyca

comm = MPI.COMM_WORLD
# end. you are free to update anything below here

# split remd.x.000 to N cores and do calc_surf in parallel
root_dir = "../../tests/data/nogit/remd/"
traj_name = root_dir + "/remd.x.000"
parm_name = root_dir + "myparm.top"

# load to TrajReadOnly
traj = io.load(traj_name, parm_name)

# mapping different chunk of `traj` in N cores
# need to provide `comm`
# save `total_arr` to rank=0
# others: total_arr = None
total_arr = pymap(comm, pyca.calc_molsurf, traj, "@CA", top=traj.top, root=0)

if comm.rank != 0:
    assert total_arr is None

if comm.rank == 0:
    # skip final array since its shape might be different from the rest
    t0 = np.asarray(total_arr[:-1]).flatten()
    t1 = np.asarray(total_arr[-1]).flatten()
    t = np.append(t0, t1)
    print ('total array len: ', t.shape[0])

    # assert to serial values
    #t2 = pyca.calc_molsurf(traj, "@CA", dtype='ndarray')
    #assert t.shape == t2.shape
    #assert np.any(t == t2) == True
