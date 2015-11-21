# (require: mpi4py, numpy)
# mpirun -n 4 python mpi_cal_molsurf_0.py

# always add those lines to your code
import numpy as np
from mpi4py import MPI
import pytraj as pt
from pytraj.testing import aa_eq

comm = MPI.COMM_WORLD
# end. you are free to update anything below here

# split remd.x.000 to N cores and do calc_surf in parallel
root_dir = "../../tests/data/"
traj_name = root_dir + "tz2.ortho.nc"
parm_name = root_dir + "tz2.ortho.parm7"

# load to TrajectoryIterator
traj = pt.iterload(traj_name, parm_name)
# print(traj)

# mapping different chunk of `traj` in N cores
# need to provide `comm`
# save `total_arr` to rank=0
# others: total_arr = None
out_parallel = pt.pmap_mpi(pt.calc_molsurf, traj, "!:WAT", top=traj.top)

if comm.rank != 0:
    assert out_parallel is None

if comm.rank == 0:
    out_serial = pt.calc_molsurf(traj, "!:WAT", dtype='ndarray')
    aa_eq(out_serial, out_parallel['MSURF_00000'])
