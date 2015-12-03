"""calculat RMSD for 8 replica trajs using 8 cores.
Reference frame is the 1st frame of remd.x.000
(embarrassingly parallel)

System: 17443 atoms, 1000 frames, netcdf, 8 replicas (000 to 007), 200Mb/replica
mpirun -n 8 python rmsd_mpi_example.py
"""

import numpy as np
from pytraj import io
from pytraj import Frame
import pytraj.common_actions as pyca
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

root_dir = "../../tests/data/nogit/remd/"
fname = root_dir + "/remd.x.00" + str(rank) # 000, 001, 002, 003 ...
top_name = root_dir + "myparm.top"

traj = io.load(fname, top_name)
n_atoms =  traj.top.n_atoms
n_frames = traj.n_frames

if rank == 0:
    ref = traj[0]
    ref_xyz = np.asarray(ref.xyz, dtype=np.float64)
else:
    ref = None
    ref_xyz = np.empty((n_atoms, 3), dtype=np.float64)

# broadcast ref_xyz to other cores from master
comm.Bcast([ref_xyz, MPI.DOUBLE])

if rank != 0:
    # need to reconstruct ref
    ref = Frame()
    ref.append_xyz(ref_xyz)

# got segmentation fault if not making a copy of ref
# TODO : check MPI stuff
_ref = ref.copy()

def rmsd_mpi(traj, _ref):
    arr0 = pyca.calc_rmsd(traj, "@CA", top=traj.top, ref= _ref)
    return arr0

arr0 = rmsd_mpi(traj, _ref)

if rank == 0:
    data = np.empty(size * traj.n_frames, dtype=np.float64)
else:
    data = None

data = comm.gather(arr0, root=0)
if rank == 0:
    all_rmsd = np.asarray(data).flatten()
    np.savetxt("mpi_rmsd.txt", all_rmsd)

    # make sure to reproduce serial version
    # YES
    #sarr = np.empty((size, traj.n_frames))
    #REF = None
    #for i in range(size):
    #    fname = "remd.x.00" + str(i)
    #    straj = io.load(fname, traj.top)
    #    if i == 0:
    #        REF = straj[0]
    #    sarr[i] = straj.calc_rmsd("@CA", REF)
    #print(np.any(all_rmsd == sarr.flatten()))
