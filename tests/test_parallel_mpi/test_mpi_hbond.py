import pytraj as pt
from mpi4py import MPI
from pytraj.testing import aa_eq

rank = MPI.COMM_WORLD.rank

traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
hbond_data_serial = pt.search_hbonds(traj, dtype='dict')
hbond_data_pmap_mpi = pt.pmap_mpi(pt.search_hbonds, traj)

if rank == 0:
    assert sorted(hbond_data_serial.keys()) == sorted(hbond_data_pmap_mpi.keys())
    for key in hbond_data_serial.keys():
        aa_eq(hbond_data_serial[key], hbond_data_pmap_mpi[key])
