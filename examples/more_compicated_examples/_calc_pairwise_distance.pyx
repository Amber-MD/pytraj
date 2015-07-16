from cpython.array cimport array
import math
import numpy as np
cimport numpy as np
from pytraj.base import *
from pytraj.common_actions import distance

def calc_pairwise_distance():
    cdef int i, j
    cdef double[:, ::1] arr
    traj = FrameArray(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")

    # extract 11th frame (index start from 0)
    frame0 = traj[10]

    print distance(frame0.coords[0:3], frame0.coords[96:99])
    natoms = frame0.n_atoms
    arr = np.empty((natoms, natoms))

    # this is the demo. For large system, use cpptraj or Cython
    # for nested loop (which is quite expensive in Python)
    for i in range(frame0.n_atoms):
        for j in range(frame0.n_atoms):
            arr[i, j] = distance(frame0.atom_xyz(i), frame0.atom_xyz(j))

from time import time
t0 = time()
calc_pairwise_distance()
print time() - t0
