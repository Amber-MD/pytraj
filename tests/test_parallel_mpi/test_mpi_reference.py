#!/usr/bin/env python

import sys
import pytraj as pt
from pytraj.testing import aa_eq

try:
    from mpi4py import MPI
except ImportError:
    print('skip if not having mpi4py')
    sys.exit(0)

comm = MPI.COMM_WORLD

rank = comm.rank

traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
data_ref3 = pt.pmap_mpi(pt.rmsd, traj, '@CA', ref=3)
# data_ref0 = pt.pmap_mpi(pt.rmsd, traj, '@CA')

if rank == 0:
    # saved_data_ref0 = pt.rmsd(traj, '@CA')
    saved_data_ref3 = pt.rmsd(traj, '@CA', ref=3)

    # aa_eq(data_ref0['RMSD_00001'], saved_data_ref0)
    aa_eq(data_ref3['RMSD_00001'], saved_data_ref3)
