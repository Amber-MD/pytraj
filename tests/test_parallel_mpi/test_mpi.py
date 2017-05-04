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
def test_mpi():
    comm = MPI.COMM_WORLD

    rank = comm.rank

    traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
    data = pt.pmap_mpi(pt.radgyr, traj, '*')

    frame_indices = range(10, 50)
    data_frame_indices = pt.pmap_mpi(
        pt.radgyr, traj, '*', frame_indices=frame_indices)
    data_frame_indices_cpp_style = pt.pmap_mpi(
        [
            'radgyr nomax',
        ], traj, frame_indices=frame_indices)

    if rank == 0:
        saved_data = pt.radgyr(traj, '*')
        saved_data_frame_indices = pt.radgyr(
            traj, '*', frame_indices=frame_indices)
        aa_eq(data['RoG_00000'], saved_data)
        aa_eq(data_frame_indices['RoG_00000'], saved_data_frame_indices)
        aa_eq(data_frame_indices_cpp_style['RoG_00000'],
              saved_data_frame_indices)


if __name__ == '__main__':
    test_mpi()
