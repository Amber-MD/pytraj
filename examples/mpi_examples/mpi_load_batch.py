# (require: mpi4py, numpy)
# mpirun -n 4 python mpi_cal_molsurf_0.py

# always add those lines to your code
import pytraj as pt
from mpi4py import MPI
comm = MPI.COMM_WORLD
from pytraj.parallel import _load_batch_pmap

# end. you are free to update anything below here
from pytraj.testing import aa_eq

root_dir = "../../tests/data/nogit/tip3p/"
traj_name = root_dir + "md.nc"
parm_name = root_dir + "tc5bwat.top"

# load to TrajectoryIterator
traj = pt.iterload(traj_name, parm_name, frame_slice=(0, 4000))

# make a list of things you want to
lines = ['autoimage', 'distance :3 :18', 'molsurf @CA']

# gather the data to 1st core (rank=0)
n_cores = comm.size
total_arr = _load_batch_pmap(n_cores, traj, lines=lines, mode='mpi')

if comm.rank != 0:
    assert total_arr is None

if comm.rank == 0:
    from pytraj.tools import concat_dict
    print(concat_dict(x[1] for x in total_arr))
