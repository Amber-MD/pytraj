#!/usr/bin/env python

import sys
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
def test_mpi_use_reference():

    comm = MPI.COMM_WORLD

    rank = comm.rank

    traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
    data_ref3 = pt.pmap_mpi(pt.rmsd, traj, '@CA', ref=3)
    # data_ref0 = pt.pmap_mpi(pt.rmsd, traj, '@CA')

    if rank == 0:
        # saved_data_ref0 = pt.rmsd(traj, '@CA')
        saved_data_ref3 = pt.rmsd(traj, '@CA', ref=3)

        # aa_eq(data_ref0['RMSD_00001'], saved_data_ref0)
        aa_eq(data_ref3['RMSD_00001'], saved_data_ref3)


if __name__ == '__main__':
    test_mpi_use_reference()
