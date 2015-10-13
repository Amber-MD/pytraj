# require: mpi4py
# mpirun -n 4 python mpi_cal_molsurf_0.py

# always add those lines to your code
import numpy as np
from mpi4py import MPI
from pytraj.parallel import map_mpi as pymap

# load pytraj
import pytraj as pt

comm = MPI.COMM_WORLD
# end. you are free to update anything below here

# split remd.x.000 to N cores and do calc_surf in parallel
root_dir = "../../tests/data/"
traj_name = root_dir + "/tz2.ortho.nc"
parm_name = root_dir + "/tz2.ortho.parm7"

# load to TrajectoryIterator
traj = pt.iterload(traj_name, parm_name)

# mapping different chunk of `traj` in N cores
# need to provide `comm`
# save `total_arr` to rank=0
# others: total_arr = None
total_arr = pymap(pt.search_hbonds, traj, ':1-13', dtype='dict')

if comm.rank != 0:
    assert total_arr is None

if comm.rank == 0:
    # save data to disk to read later
    # read: pt.read_pickle('output/hbond_data.pk')
    pt.to_pickle(total_arr, 'output/hbond_data.pk')
