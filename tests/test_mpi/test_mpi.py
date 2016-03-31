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
data = pt.pmap_mpi(pt.radgyr, traj, '*')

frame_indices = range(10, 50)
data_frame_indices = pt.pmap_mpi(pt.radgyr,
                                 traj,
                                 '*',
                                 frame_indices=frame_indices)
data_frame_indices_cpp_style = pt.pmap_mpi(
    ['radgyr nomax', ],
    traj,
    frame_indices=frame_indices)

if rank == 0:
    saved_data = pt.radgyr(traj, '*')
    saved_data_frame_indices = pt.radgyr(traj,
                                         '*',
                                         frame_indices=frame_indices)
    aa_eq(data['RoG_00000'], saved_data)
    aa_eq(data_frame_indices['RoG_00000'], saved_data_frame_indices)
    aa_eq(data_frame_indices_cpp_style['RoG_00000'], saved_data_frame_indices)
