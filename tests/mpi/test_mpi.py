#!/usr/bin/env python

import pytraj as pt
from pytraj.testing import aa_eq
from mpi4py import MPI

comm = MPI.COMM_WORLD

rank = comm.rank

traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
data = pt.pmap_mpi(pt.radgyr, traj, '*')

if rank == 0:
    data = pt.tools.flatten(data)
    saved_data = pt.radgyr(traj, '*')
    aa_eq(data, saved_data)
