# (require: mpi4py, numpy)
# mpirun -n 4 python mpi_cal_molsurf_0.py

# always add those lines to your code
from mpi4py import MPI
from pytraj.parallel import map as pymap
from pytraj import io
import pytraj.common_actions as pyca

comm = MPI.COMM_WORLD
# end. you are free to update anything below here

# split remd.x.000 to N cores and do calc_surf in parallel
root_dir = "../tests/data/nogit/remd/"
traj_name = root_dir + "/remd.x.000"
parm_name = root_dir + "myparm.top"

# load to TrajReadOnly
traj = io.load(traj_name, parm_name)

# mapping different chunk of `fa` in N cores
# need to provide `comm`
arr = pymap(comm, fa, pyca.calc_molsurf, "@CA")
print ("rank = %s, return arr with len=%s" % (comm.rank, len(arr)))
