import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq

try:
    from mpi4py import MPI
    has_mpi4py = True
except ImportError:
    MPI = None
    has_mpi4py = False


@unittest.skipUnless(has_mpi4py, 'must have mpi4py')
def test_mpi_hbond():
    rank = MPI.COMM_WORLD.rank

    traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
    hbond_data_serial = pt.search_hbonds(traj, dtype='dict')
    hbond_data_pmap_mpi = pt.pmap_mpi(pt.search_hbonds, traj)

    if rank == 0:
        assert sorted(hbond_data_serial.keys()) == sorted(
            hbond_data_pmap_mpi.keys())
        for key in hbond_data_serial.keys():
            aa_eq(hbond_data_serial[key], hbond_data_pmap_mpi[key])


if __name__ == '__main__':
    test_mpi_hbond()
