#!/usr/bin/env python

import pytraj as pt
from pytraj.testing import aa_eq
from mpi4py import MPI

comm = MPI.COMM_WORLD

rank = comm.rank

traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
data = pt.pmap_mpi(pt.radgyr, traj, '*')

if rank == 0:
    saved_data = pt.radgyr(traj, '*')
    aa_eq(data['RoG_00000'], saved_data)
