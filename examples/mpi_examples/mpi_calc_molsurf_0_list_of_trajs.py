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
root_dir = "../tests/data/nogit/remd/"
parm_name = root_dir + "myparm.top"
top = io.load(parm_name)

# load to TrajReadOnly
trajlist = []
for i in range(comm.size):
    ext = "00" + str(i)
    traj_name = root_dir + "/remd.x." + ext # 000, 001, 002
    trajlist.append(io.load(traj_name, top))

# mapping different traj to N cores
# need to provide `comm`
print (len(trajlist))
arr = pymap(comm, trajlist, pyca.calc_molsurf, "@CA", top=top)
print ("rank = %s, return arr with len=%s" % (comm.rank, len(arr)))

# gathering the data to root=0
#if comm.rank == 0:
#    total_arr =  np.empty(comm.size)
#else:
#    total_arr = None
total_arr = comm.gather(arr, root=0)

if comm.rank == 0:
    # skip final array since its shape might be different from the rest
    t0 = np.asarray(total_arr[:-1]).flatten()
    t1 = np.asarray(total_arr[-1]).flatten()
    t = np.append(t0, t1)
    print ('total array len: ', t.shape[0])

    # assert to serial values
    #t2 = pyca.calc_molsurf(trajlist, "@CA", dtype='ndarray')
    #assert t.shape == t2.shape
    #assert np.any(t == t2) == True
