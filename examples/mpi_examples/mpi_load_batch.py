# example:
# mpirun -n 4 your_script.py

# always add those lines to your code
import pytraj as pt
from mpi4py import MPI

# create ``comm`` so you can have the info about n_cpus, cpu id
comm = MPI.COMM_WORLD

# you are free to update anything below here
# _load_batch_pmap is temp method name, will be changed in future
from pytraj.parallel import _load_batch_pmap

root_dir = "../../tests/data/nogit/tip3p/"
traj_name = root_dir + "md.nc"
parm_name = root_dir + "tc5bwat.top"

# load to TrajectoryIterator
traj = pt.iterload(traj_name, parm_name, frame_slice=(0, 4000))

# make a list of things you want to
lines = ['autoimage', 'distance :3 :18', 'molsurf @CA']

# gather the data to 1st core (rank=0)
# 
n_cores = comm.size
data = _load_batch_pmap(n_cores, traj, lines=lines, mode='mpi', dtype='dict')

if comm.rank != 0:
    assert data is None

if comm.rank == 0:
    from pytraj.tools import concat_dict
    # each core return a tuple (core_id, dict)
    # so you need to concat the dict
    # use `from pytraj.tools import concat_dict
    print(concat_dict(x[1] for x in data))
